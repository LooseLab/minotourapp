# Program to parse fastq files generated by ont2d
import json
import os
import platform  # MS
import queue
import sys
import threading
import time

# import MySQLdb #note problem installing on python 3
import configargparse
import dateutil.parser
import numpy as np
import requests
from Bio import SeqIO
from watchdog.events import FileSystemEventHandler
from watchdog.observers.polling import PollingObserver as Observer
from watchdog.utils.compat import queue


class Urls(object):
    """ Holds the list of api endpoints.
    """

    def __init__(self, hostname, port):
        self.baseurl = 'http://{}:{}/api/v1'.format(hostname, port)

        self.RUNS = self.baseurl + '/runs/'
        self.READ_TYPES = self.baseurl + '/readtypes/'


class FastqDetails(object):
    """ Holds the header information about a fastq read.
    """

    def __init__(self, start_time, barcode, runid, read, channel):
        self.start_time = start_time
        self.runid = runid
        self.read = read
        self.channel = channel

        if barcode:
            self.is_barcoded = True
            self.barcode = barcode
        else:
            self.is_barcoded = False
            self.barcode = 'unclassified'

    def __str__(self):
        message = ("FastqDetails [start_time: {}, is_barcoded: {}, runid: {}, "
                   "read: {}, channel: {}]")

        return message.format(
            self.start_time,
            self.is_barcoded,
            self.runid,
            self.read,
            self.channel)


def parse_fastq_file(fastq, rundict):

    print('Parsing file {}'.format(fastq))

    for record in SeqIO.parse(fastq, "fastq"):

        fastq_details = parse_fastq_read_header(record.description)

        if fastq_details.runid not in rundict:
            rundict.update({
                fastq_details.runid: Runcollection(ARGS)
            })

            rundict[fastq_details.runid].add_run(fastq_details)

        rundict[fastq_details.runid].add_read(record, fastq_details)


def parse_fastq_read_header(description):
    descriptiondict = dict()
    descriptors = description.split(" ")
    del descriptors[0]
    for item in descriptors:
        bits = item.split("=")
        descriptiondict[bits[0]] = bits[1]

    fastq_details = FastqDetails(
        descriptiondict['start_time'],
        descriptiondict['barcode'],
        descriptiondict['runid'],
        descriptiondict['read'],
        descriptiondict['ch']
    )

    return fastq_details


class Runcollection():
    """ Holds information about a particular run.
    """

    def __init__(self, args):
        self.args = args
        self.readid = dict()
        self.readnames = list()
        self.cumulength = 0
        self.readcount = 0
        self.readlengths = list()
        self.timeid = dict()
        self.chandict = list()
        self.runidlink = ""
        self.readtypes = dict()
        self.statsrecord = dict()
        self.barcodes = dict()

        # update list of read types
        readtypes = requests.get(URLS.READ_TYPES, headers=header)

        for readtype in json.loads(readtypes.text):
            name = readtype["name"]
            url = readtype["url"]

            self.readtypes.update({
                name: url
            })

    def add_run(self, fastq_details):
        print("Seen a new run")
        # Test to see if the run exists
        r = requests.get(URLS.RUNS, headers=header)
        print(r.text)
        print(type(r.text))

        if fastq_details.runid not in r.text:
            print('Need to create this run.')
            #We need to come up with a way of identifying the run name.
            runname = self.args.run_name

            #if "barcode" in descriptiondict.keys():
            #    is_barcoded = True
            #    barcoded = "barcoded"
            #else:
            #    is_barcoded = False
            #    barcoded = "unclassified"

            #createrun = requests.post(args.full_host+'api/v1/runs/',
            #headers=header, json={"run_name": runname, "run_id": runid,
            #"barcode": barcoded, "is_barcoded":is_barcoded})

            createrun = requests.post(
                URLS.RUNS,
                headers=header,
                json={
                    "run_name": runname,
                    "run_id": fastq_details.runid,
                    "is_barcoded": fastq_details.is_barcoded
                }
            )

            print(createrun.text)

            if createrun.status_code != 201:
                print(createrun.status_code)
                print(createrun.text)
                print("Houston - we have a problem!")
            else:
                print(json.loads(createrun.text)["url"])
                self.runidlink = json.loads(createrun.text)["url"]
        else:
            print('Run Exists.')
            for run in json.loads(r.text):
                if run["run_id"] == fastq_details.runid:
                    self.runidlink = run["url"]

            # Now fetch a list of reads that already exist at that location
            checkreads = requests.get(self.runidlink + 'readnames', headers=header)
            for readid in json.loads(checkreads.text):
                self.readnames.append(readid)

            response = requests.get(self.runidlink + 'barcodes', headers=header)
            for item in json.loads(response.text):
                self.barcodes.update({
                    item.name: {
                        'url': item.url,
                        'id': item.id,
                        'run': item.run
                    }
                })


    def add_read_db(self, readid, read, channel, barcode, sequence,
                    quality, ispass, typelink, starttime):
        #print(runid,readid,read,channel,barcode,sequence,quality,ispass,type,starttime)
        #runlink = args.full_host+'api/v1/runs/' + str(self.runidlink) + "/"
        runlink = self.runidlink
        runlinkaddread = self.runidlink + "reads/"
        #typelink = args.full_host+'api/v1/readtypes/' + str(type) + "/"
        print(runlink, typelink)

        if readid not in self.readnames:
            payload = {
                'run_id': runlink,
                'read_id': readid,
                'read': read,
                "channel": channel,
                'barcode': barcode,
                'sequence': sequence,
                'quality': quality,
                'is_pass': ispass,
                'start_time': starttime,
                'type':typelink
            }

            requests.post(runlinkaddread, headers=header, json=payload)

        else:
            print("Read Seen")

    def add_or_update_stats(self, sample_time, total_length, max_length,
                            min_length, average_length, number_of_reads,
                            number_of_channels, read_type_id):
        runlink = ARGS.full_host+'api/v1/runs/' + str(self.runidlink) + "/"
        typelink = ARGS.full_host+'api/v1/readtypes/' + str(read_type_id) + "/"
        payload = {
            'run_id': runlink,
            'sample_time': sample_time,
            'total_length': total_length,
            'max_length': max_length,
            'min_length': min_length,
            'average_length': average_length,
            'number_of_reads': number_of_reads,
            'number_of_channels': number_of_channels,
            'type': typelink
        }

        if sample_time not in self.statsrecord.keys():
            createstats = requests.post(
                ARGS.full_host+'api/v1/statistics/',
                headers=header,
                json=payload
            )

            print(json.loads(createstats.text)["id"])
            self.statsrecord[sample_time] = dict()
            self.statsrecord[sample_time]["id"] = json.loads(createstats.text)["id"]
            self.statsrecord[sample_time]['payload'] = payload
        else:
            print("Seen these stats before!")
            if self.statsrecord[sample_time]['payload'] != payload:
                print("Record changed - needs updating")
                print(payload)
                # payload['id']=self.statsrecord[sample_time]['id']
                # print(payload)
                url = "{}api/v1/statistics/{}/".format(
                    ARGS.full_host,
                    str(self.statsrecord[sample_time]['id'])
                )

                requests.delete(
                    url,
                    headers=header
                )

                createstats = requests.post(
                    ARGS.full_host + 'api/v1/statistics/',
                    headers=header,
                    json=payload
                )

                self.statsrecord[sample_time]["id"] = json.loads(createstats.text)["id"]
                self.statsrecord[sample_time]['payload'] = payload
                print(createstats.text)
            else:
                print("Record not changed.")

    def add_read(self, record, fastq_details):
        if record.id not in self.readid:
            self.readid[record.id] = fastq_details

            tm = dateutil.parser.parse(self.readid[record.id].start_time)

            print(fastq_details.start_time)

            print(tm)

            #tm = tm - datetime.timedelta(
            #    minutes=(tm.minute % 1) - 1,
            #    seconds=tm.second,
            #    microseconds=tm.microsecond
            #)

            tm = tm.replace(second=0, microsecond=0)

            print(tm)

            if tm not in self.timeid.keys():
                self.timeid[tm] = dict()
                self.timeid[tm]["cumulength"] = 0
                self.timeid[tm]["count"] = 0
                self.timeid[tm]["readlengths"] = list()
                self.timeid[tm]["chandict"] = list()

            # This is illustrating how to access the sequence but will be a memory problem
            # self.readid[record.id]["seq"]=record.seq
            # self.readid[record.id]["qual"]=record.format("qual")
            if self.readid[record.id].channel not in self.chandict:
                self.chandict.append(self.readid[record.id].channel)

            if self.readid[record.id].channel not in self.timeid[tm]["chandict"]:
                self.timeid[tm]["chandict"].append(self.readid[record.id].channel)

            if self.readid[record.id].is_barcoded:
                barcode_name = self.readid[record.id].barcode

                if barcode_name in self.barcodes.keys():
                    barcode_url = self.barcodes[barcode_name].url

                else:
                    body = {
                        'name': barcode_name,
                        'run': self.runidlink
                    }

                    response = requests.post(
                        self.runidlink + 'barcodes/',
                        headers=header,
                        json=body
                    )

                    if response.status_code != 201:
                        print(response.status_code)
                        print(response.text)
                        print("Houston - we have a problem2!")
                    else:
                        print(json.loads(response.text)["url"])
                        barcode_url = json.loads(response.text)["url"]

            else:
                barcode_url = None

            self.add_read_db(
                record.id,
                self.readid[record.id]["read"],
                self.readid[record.id]["ch"],
                barcode_url,
                str(record.seq),
                record.format('fastq').split('\n')[3],
                True,
                self.readtypes["Template"],
                self.readid[record.id]["start_time"]
            )

            self.readid[record.id]["len"] = len(record.seq)
            self.cumulength += len(record.seq)
            self.timeid[tm]["cumulength"] += len(record.seq)
            self.readcount += 1
            self.timeid[tm]["count"] += 1
            self.readlengths.append(len(record.seq))
            self.timeid[tm]["readlengths"].append(len(record.seq))
            # print(record.id)
            # print(record.seq)
            # print(record.description)
            # print(record.format("qual"))
            # print(len(record.seq))

    def read_count(self):
        return len(self.readid)

    def mean_median_std_max_min(self):
        return np.average(self.readlengths),\
            np.median(self.readlengths),\
            np.std(self.readlengths),\
            np.max(self.readlengths),\
            np.min(self.readlengths)

    def parse1minwin(self):
        for my_time in sorted(self.timeid):
            self.add_or_update_stats(
                str(my_time),
                self.timeid[my_time]["cumulength"],
                np.max(self.timeid[my_time]["readlengths"]),
                np.min(self.timeid[my_time]["readlengths"]),
                np.around(np.average(self.timeid[my_time]["readlengths"]), decimals=2),
                self.timeid[my_time]["count"],
                len(self.timeid[my_time]["chandict"]),
                self.readtypes["Template"]
            )
#print(time, self.timeid[time]["cumulength"], np.max(self.timeid[time]["readlengths"]),
#      np.min(self.timeid[time]["readlengths"]), np.average(self.timeid[time]["readlengths"]),
#      self.timeid[time]["count"], len(self.timeid[time]["chandict"]))


class MyHandler(FileSystemEventHandler):
    """ To come.
    """

    num_threads = 10

    def __init__(self):
        """ Collect information about files already in the folders.
        """

        self.fastq_files_queue = queue.Queue(maxsize=0)

        #self.args = args
        #self.creates = get_existing_fastq_files(args.watchdir)
        self.creates = dict()
        self.rundict = dict()

        self.process_queue()

        #t = threading.Thread(target=self.processfiles)

        #try:
        #    t.start()
        #    print("Watchdog started")

        #except (KeyboardInterrupt, SystemExit):
        #    t.stop()

    def process_queue(self):
        for i in range(self.num_threads):
            worker = threading.Thread(target=self.process_queue_items)
            worker.setDaemon(True)
            worker.start()

        self.fastq_files_queue.join()

    def process_queue_items(self):
        while True:
            item = self.fastq_files_queue.get()
            print(item)
            self.processfiles(item[0])
            self.fastq_files_queue.task_done()

    def processfiles(self, filename):
        """ To come.
        """

        parse_fastq_file(filename, self.rundict)

    def processfiles2(self):
        """ To come.
        """

        while True:
            print("processfiles running")

            for fastqfile, createtime in sorted(self.creates.items(), key=lambda x: x[1]):
                delaytime = 0

                """ file created 5 sec ago, so should be complete.
                    For simulations we make the time longer.
                """
                if (int(createtime) + delaytime < time.time()):
                    print(fastqfile)
                    del self.creates[fastqfile]
                    parse_fastq_file(fastqfile, self.rundict)

            time.sleep(5)

    def process_event(self, event):
        if event.src_path.endswith('fastq'):
            print('event occured {}'.format(event))
            self.fastq_files_queue.put((event.src_path, time.time()))
            self.creates[event.src_path] = time.time()

    def on_created(self, event):
        self.process_event(event)

    def on_modified(self, event):
        self.process_event(event)

def get_existing_fastq_files(path, handler):
    """ Retrives all fastq files that already exist in the watch dir.
        Watchdog will only monitor for new fastq files.
    """

    file_list_dict = dict()

    if os.path.isdir(path):

        #print("caching existing fastq files in: %s" % (path))

        for path, dirs, files in os.walk(path):

            for the_file in files:

                if (the_file.endswith(".fastq")):
                    print(the_file)
                    #file_list_dict[os.path.join(path, the_file)] = os.stat(os.path.join(path, the_file)).st_mtime
                    handler.fastq_files_queue.put((
                        os.path.join(path, the_file),
                        os.stat(os.path.join(path, the_file)).st_mtime
                    ))

    print("Found {} existing fast5 files to process first.".format(handler.fastq_files_queue.qsize()))


# noinspection PyGlobalUndefined
if __name__ == '__main__':

    OPER = platform.system()
    if OPER is 'Windows':  # MS
        OPER = 'windows'
    else:
        OPER = 'linux'  # MS

    if OPER is 'linux':
        config_file = os.path.join(os.path.sep, \
            os.path.dirname(os.path.realpath('__file__' \
            )), 'minfq_posix.config')

    if OPER is 'windows':
        config_file = os.path.join(os.path.sep, sys.prefix,
            'minfq_windows.config')

    parser = \
        configargparse.ArgParser(description='minFQ: A program to analyse minION fastq files in real-time or post-run.' \
                                 , default_config_files=[config_file])

    parser.add(
        '-w',
        '--watch-dir',
        type=str,
        required=True,
        default=None,
        help='The path to the folder containing the downloads directory with fast5 reads to analyse - e.g. C:\\data\\minion\\downloads (for windows).',
        dest='watchdir'
    )

    parser.add(
        '-n',
        '--name',
        type=str,
        required=True,
        default=None,
        help='The run name you wish to provide.',
        dest='run_name',
    )

    parser.add(
        '-k',
        '--key',
        type=str,
        required=True,
        default=None,
        help='The api key for uploading data.',
        dest='api_key',
    )

    parser.add(
        '-p',
        '--port',
        type=int,
        #required=True,
        default=8100,
        help='The port number for the local server.',
        dest='port_number',
    )

    parser.add(
        '-hn',
        '--hostname',
        type=str,
        #required=True,
        default='127.0.0.1',
        help='The run name you wish to provide.',
        dest='host_name',
    )

    ARGS = parser.parse_args()

    ARGS.full_host = "http://" + ARGS.host_name + ":" + str(ARGS.port_number) + "/"

    URLS = Urls(ARGS.host_name, ARGS.port_number)

    #GLobal creation of header (needs fixing)
    header = {
        'Authorization': 'Token ' + ARGS.api_key,
        'Content-Type': 'application/json'
    }

    print(ARGS.watchdir)

    event_handler = MyHandler()
    observer = Observer()
    observer.schedule(event_handler, path=ARGS.watchdir, recursive=True)
    observer.daemon = True

    observer.start()

    get_existing_fastq_files(ARGS.watchdir, event_handler)

    try:
        while True:
            time.sleep(1)

    except (KeyboardInterrupt, SystemExit):
        print("catching a ctrl-c event")
        # my_client.stop()
        observer.stop()
        observer.join()
        # die_nicely(oper)
