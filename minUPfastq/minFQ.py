import datetime
import json
import os
import platform  # MS
import sys
import threading
import time

# import MySQLdb #note problem installing on python 3
import configargparse
import dateutil.parser
import numpy as np
import requests
from Bio import SeqIO
from tqdm import tqdm
from watchdog.events import FileSystemEventHandler
from watchdog.observers.polling import PollingObserver as Observer


def parsefastq(fastq, rundict):
    #print ('processing reads')
    for record in SeqIO.parse(fastq, "fastq"):
        descriptiondict = parsedescription(record.description)
        if descriptiondict["runid"] not in rundict:
            rundict[descriptiondict["runid"]] = Runcollection(args)
            rundict[descriptiondict["runid"]].add_run(descriptiondict)
        rundict[descriptiondict["runid"]].add_read(record, descriptiondict,fastq)


def parsedescription(description):
    descriptiondict = dict()
    descriptors = description.split(" ")
    del descriptors[0]
    for item in descriptors:
        bits = item.split("=")
        descriptiondict[bits[0]] = bits[1]
    return descriptiondict


class Runcollection():

    def __init__(self, args):
        self.args = args
        self.readid = dict()
        self.readnames=list()
        self.cumulength = 0
        self.readcount = 0
        self.readlengths = list()
        self.timeid = dict()
        self.chandict = list()
        self.runidlink = ""
        self.readtypes=dict()
        self.statsrecord=dict()
        self.barcodes = dict()
        self.runid = 0

    def get_readnames_by_run(self):
        content = requests.get(self.runidlink + 'readnames', headers=header)

        content_json = json.loads(content.text)

        number_pages = content_json.get('number_pages')

        print ('Requesting readnames already uploaded to databases.')

        for page in tqdm(range(number_pages)):

            url = self.runidlink + 'readnames?page=%s'.format(page)

            content = requests.get(url, headers=header)

            for read in json.loads(content.text):
                self.readnames.append(read)

    def add_run(self, descriptiondict):
        print("Seen a new run")
        # Test to see if the run exists
        r = requests.get(args.full_host+'api/v1/runs', headers=header)

        runid = descriptiondict["runid"]
        if runid not in r.text:
            #print ('Need to create this run.')
            #We need to come up with a way of identifying the run name.
            runname = self.args.run_name
            if "barcode" in descriptiondict.keys():
                is_barcoded = True
                barcoded = "barcoded"
            else:
                is_barcoded = False
                barcoded = "unclassified"
            createrun = requests.post(args.full_host+'api/v1/runs/', headers=header, json={"run_name": runname, "run_id": runid, "barcode": barcoded, "is_barcoded":is_barcoded})

            if createrun.status_code != 201:
                print (createrun.status_code)
                print (createrun.text)
                print ("Houston - we have a problem!")
            else:
                self.runidlink = json.loads(createrun.text)["url"]
                self.runid = json.loads(createrun.text)["id"]
        else:
            #print ('Run Exists.')
            for run in json.loads(r.text):
                if run["run_id"] == runid:
                    print (run)
                    self.runidlink = run["url"]
                    self.runid = run["id"]

            #
            # Now fetch a list of reads that already exist at that location
            #
            self.get_readnames_by_run()

        readtypes=requests.get(args.full_host+'api/v1/readtypes', headers=header)

        for readtype in json.loads(readtypes.text):
            self.readtypes[readtype["name"]]=readtype["url"]

        response = requests.get(
            str(self.runidlink) + "barcodes/",
            headers=header
        )

        for item in json.loads(response.text):
            self.barcodes.update({
                item['name']: item['url']
            })

    def add_read_db(self,runid,readid,read,channel,barcode,sequence,quality,ispass,type,starttime):
        #print(runid,readid,read,channel,barcode,sequence,quality,ispass,type,starttime)
        #runlink = args.full_host+'api/v1/runs/' + str(self.runidlink) + "/"
        runlink = self.runidlink
        runlinkaddread = self.runidlink + "reads/"
        #typelink = args.full_host+'api/v1/readtypes/' + str(type) + "/"
        typelink = type
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
                'type': typelink
            }
            createread = requests.post(runlinkaddread, headers=header, json=payload)
        else:
            print ("Read Seen")

    def add_or_update_stats(self,sample_time,total_length,max_length,min_length,average_length,number_of_reads,number_of_channels,type):
        runlink = args.full_host+'api/v1/runs/' + str(self.runidlink) + "/"
        typelink = args.full_host+'api/v1/readtypes/' + str(type) + "/"
        payload = {'run_id': runlink, 'sample_time':sample_time,'total_length':total_length,'max_length':max_length,'min_length':min_length,'average_length':average_length,'number_of_reads':number_of_reads,'number_of_channels':number_of_channels,'type':typelink}
        if sample_time not in self.statsrecord.keys():
            createstats=requests.post(args.full_host+'api/v1/statistics/', headers=header, json=payload)
            print (json.loads(createstats.text)["id"])
            self.statsrecord[sample_time]=dict()
            self.statsrecord[sample_time]["id"]=json.loads(createstats.text)["id"]
            self.statsrecord[sample_time]['payload']=payload
        else:
            print ("Seen these stats before!")
            if self.statsrecord[sample_time]['payload'] != payload:
                print ("Record changed - needs updating")
                print(payload)
                # payload['id']=self.statsrecord[sample_time]['id']
                # print (payload)
                deletestats = requests.delete(
                    args.full_host+'api/v1/statistics/' + str(self.statsrecord[sample_time]['id']) + '/',
                    headers=header)
                createstats = requests.post(args.full_host+'api/v1/statistics/', headers=header, json=payload)
                self.statsrecord[sample_time]["id"] = json.loads(createstats.text)["id"]
                self.statsrecord[sample_time]['payload'] = payload
                print(createstats.text)
            else:
                print ("Record not changed.")

    def update_read_type(self,read_id,type):
        payload = {'type': type}
        updateread = requests.patch(args.full_host+'api/v1/runs/' + str(self.runid) + "/reads/" + str(read_id) +'/', headers=header, json=payload)
        print (updateread.text)

    def check_1d2(self,readid):
        if len(readid) > 64:
            return True

    def check_pass(self,path):
        folders = os.path.split(path)
        #print folders[0]
        if 'pass' in folders[0]:
            return True
        elif 'fail' in folders[0]:
            return False
        else:
            return True #This assumes we have been unable to find either pass or fail and thus we assume the run is a pass run.

    def add_read(self, record, descriptiondict,fastq):
        passstatus=(self.check_pass(fastq))

        if record.id not in self.readid:
            self.readid[record.id] = dict()
            for item in descriptiondict:
                self.readid[record.id][item] = descriptiondict[item]

            # print self.readid[record.id]
            tm = dateutil.parser.parse(self.readid[record.id]["start_time"])
            # print tm
            tm = tm - datetime.timedelta(minutes=(tm.minute % 1) - 1,
                                         seconds=tm.second,
                                         microseconds=tm.microsecond)
            # print tm
            if tm not in self.timeid.keys():
                self.timeid[tm] = dict()
                self.timeid[tm]["cumulength"] = 0
                self.timeid[tm]["count"] = 0
                self.timeid[tm]["readlengths"] = list()
                self.timeid[tm]["chandict"] = list()

            # This is illustrating how to access the sequence but will be a memory problem
            # self.readid[record.id]["seq"]=record.seq
            # self.readid[record.id]["qual"]=record.format("qual")
            if self.readid[record.id]["ch"] not in self.chandict:
                self.chandict.append(self.readid[record.id]["ch"])

            if self.readid[record.id]["ch"] not in self.timeid[tm]["chandict"]:
                self.timeid[tm]["chandict"].append(self.readid[record.id]["ch"])

            if "barcode" in self.readid[record.id].keys() or args.cust_barc=='oddeven':

                if "barcode" in self.readid[record.id]:
                    barcode_local = self.readid[record.id]["barcode"]
                else:
                    barcode_local = ""

                if args.cust_barc=='oddeven':
                    if int(self.readid[record.id]["ch"]) % 4 == 0:
                        barcode_local = barcode_local + 'even'
                    else:
                        barcode_local = barcode_local + 'odd'

                if barcode_local not in self.barcodes.keys():

                    print(">> Found new barcode {} for run {}.".format(barcode_local, self.runidlink))

                    request_body = {
                        'name': barcode_local,
                        'run': str(self.runidlink)
                    }

                    response = requests.post(
                        str(self.runidlink) + "barcodes/",
                        headers=header,
                        json=request_body
                    )

                    if response.status_code == 201:
                        item = json.loads(response.text)

                        self.barcodes.update({
                            item['name']: item['url']
                        })

                        print(">> Barcode {} for run {} created with success.".format(item['url'], self.runidlink))

                barcode_url = self.barcodes[barcode_local]

            else:
                barcode_url = self.barcodes["No barcode"]

            if self.check_1d2(record.id):
                print ("Seen a 1D^2 read. Exiting.")
                print (record.id)
                print (len(record.id))
                firstread, secondread = record.id[:len(record.id) // 2], record.id[len(record.id) // 2:]
                print (firstread,secondread)
                self.update_read_type(secondread,self.readtypes["Complement"])
                #So here we need to a) add the new 2D read as a 2D read - then update the read status of the second read.
                self.add_read_db(
                    self.runidlink,
                    record.id,
                    self.readid[record.id]["read"],
                    self.readid[record.id]["ch"],
                    barcode_url,
                    str(record.seq),
                    record.format('fastq').split('\n')[3],
                    passstatus,
                    self.readtypes["1D^2"],
                    self.readid[record.id]["start_time"]
                )

            else:
                self.add_read_db(
                    self.runidlink,
                    record.id,
                    self.readid[record.id]["read"],
                    self.readid[record.id]["ch"],
                    barcode_url,
                    str(record.seq),
                    record.format('fastq').split('\n')[3],
                    passstatus,
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
        return np.average(self.readlengths), np.median(self.readlengths), np.std(self.readlengths), np.max(
            self.readlengths), np.min(self.readlengths)

    def parse1minwin(self):
        for time in sorted(self.timeid):
            self.add_or_update_stats(str(time),self.timeid[time]["cumulength"], np.max(self.timeid[time]["readlengths"]),
                  np.min(self.timeid[time]["readlengths"]), np.around(np.average(self.timeid[time]["readlengths"]), decimals=2),
                  self.timeid[time]["count"], len(self.timeid[time]["chandict"]),self.readtypes["Template"])
            #print(time, self.timeid[time]["cumulength"], np.max(self.timeid[time]["readlengths"]),
            #      np.min(self.timeid[time]["readlengths"]), np.average(self.timeid[time]["readlengths"]),
            #      self.timeid[time]["count"], len(self.timeid[time]["chandict"]))


def file_dict_of_folder_simple(path):
    file_list_dict = dict()
    # ref_list_dict=dict()
    print("File Dict Of Folder Called")
    counter = 0
    if os.path.isdir(path):
        print("caching existing fastq files in: %s" % (path))
        for path, dirs, files in tqdm(os.walk(path)):
            for f in files:
                # print f
                counter += 1
                #print(counter)
                # if (("downloads" in path )):
                # if ("muxscan" not in f and args.callingdir not in path and f.endswith(".fast5") ):
                if (f.endswith(".fastq")):
                    file_list_dict[os.path.join(path, f)] = os.stat(os.path.join(path, f)).st_mtime
                    # try:
                    #    file_descriptor = update_file_descriptor(os.path.join(path, f),file_descriptor)
                    # except:
                    #    pass
    print("processed %s files" % (counter))
    print("found %d existing fastq files to process first." % (len(file_list_dict)))
    return file_list_dict


class MyHandler(FileSystemEventHandler):
    def __init__(self, args):
        """Collect information about files already in the folders"""
        self.file_descriptor = dict()
        self.args = args
        # adding files to the file_descriptor is really slow - therefore lets skip that and only update the files when we want to basecall thread_number
        # self.creates,self.file_descriptor=file_dict_of_folder(args.watchdir, self.file_descriptor)
        self.creates = file_dict_of_folder_simple(args.watchdir)
        self.processing = dict()
        self.running = True
        self.rundict = dict()
        t = threading.Thread(target=self.processfiles)

        try:
            t.start()
            print("Watchdog started")
        #except (KeyboardInterrupt, SystemExit):
        except KeyboardInterrupt:
            print ("Seen a ctrl-c")
            t.stop()
            raise

    def lencreates(self):
        return len(self.creates)

    def lenprocessed(self):
        return len(self.processed)

    def processfiles(self):
        everyten = 0
        while self.running:
            print("processfiles running")
            for fastqfile, createtime in tqdm(sorted(self.creates.items(), key=lambda x: x[1])):
                delaytime = 0
                if (int(
                        createtime) + delaytime < time.time()):  # file created 5 sec ago, so should be complete. For simulations we make the time longer.
                    #print(fastqfile)
                    del self.creates[fastqfile]
                    parsefastq(fastqfile, self.rundict)

            for runid in self.rundict:
                print("RunID", runid)
                print("Read Number:", self.rundict[runid].readcount, "Total Length:", self.rundict[runid].cumulength,
                      "Average Length", self.rundict[runid].cumulength / self.rundict[runid].readcount, "Chan Count",
                      len(self.rundict[runid].chandict))
                (mean, median, std, maxval, minval) = self.rundict[runid].mean_median_std_max_min()
                print("mean", mean, "median", median, "std", std, "max", maxval, "min", minval)
                # print self.rundict[runid].timeid
                #self.rundict[runid].parse1minwin()

            time.sleep(5)

    def on_created(self, event):
        """Watchdog counts a new file in a folder it is watching as a new file"""
        """This will add a file which is added to the watchfolder to the creates and the info file."""
        if (event.src_path.endswith(".fastq")):
            # print "seen a file", event.src_path
            self.creates[event.src_path] = time.time()
            # elif ("albacore" in event.src_path and args.callingdir not in event.src_path and args.finishdir not in event.src_path and event.src_path.endswith(".fast5")):
            #    self.albacorefiles[event.src_path] = time.time()
            #    print "seen an albacore file created"

            # self.total[event.src_path] = time.time()

    def on_modified(self, event):
        if (event.src_path.endswith(".fastq")):
            # print "seen a file", event.src_path
            self.creates[event.src_path] = time.time()

    def on_moved(self, event):
        """Watchdog considers a file which is moved within its domain to be a move"""
        """When a file is moved, we just want to update its location in the master dictionary."""
        # print "On Moved Called"
        if (event.dest_path.endswith(".fastq")):
            print("seen a fastq file move")
            # print getfilename(event.src_path), event.src_path,event.dest_path
            # try:
            # self.creates[event.dest_path] = self.creates[event.src_path]
            # del self.creates[event.src_path]
            # print "file finished basecalling:", getfilename(event.dest_path)
            #    self.albacoredone[event.dest_path] = time.time()
            # del self.albacorefiles[event.src_path]
            # del self.creates[event.src_path]
            # except Exception as e:
            #    print "Error deleting file:",e
            #    pass
            # else:
            # print "seen a file move"
            # print getfilename(event.src_path), event.src_path,event.dest_path
            # print self.creates
            # del self.creates[]

    def on_deleted(self, event):
        print("On Deleted Called", event.src_path)
        # if event.src_path in self.creates.keys():
        #     del self.creates[event.src_path]
        # """ We need to clean out any other references to these files"""
        # if getfilename(event.src_path) in self.processing.keys():
        #     del self.processing[getfilename(event.src_path)]
        # if getfilename(event.src_path) in self.processed.keys():
        #     del self.processed[getfilename(event.src_path)]


# noinspection PyGlobalUndefined
if __name__ == '__main__':

    global OPER

    OPER = platform.system()
    if OPER is 'Windows':  # MS
        OPER = 'windows'
    else:
        OPER = 'linux'  # MS
    print(OPER)  # MS

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

    parser.add(
        '-c',
        '--custom_barcode',
        type=str,
        required=False,
        default=None,
        help='Optionally split reads based on odd/even channel description. Not a standard option.',
        dest='cust_barc'
    )

    args = parser.parse_args()

    args.full_host = "http://" + args.host_name + ":" + str(args.port_number) + "/"

    #GLobal creation of header (needs fixing)
    global header
    header = {'Authorization': 'Token ' + args.api_key, 'Content-Type': 'application/json'}

    event_handler = MyHandler(args)
    observer = Observer()
    observer.schedule(event_handler, path=args.watchdir, recursive=True)
    observer.daemon = True
    observer.start()

    try:
        while 1:
            time.sleep(1)
    except KeyboardInterrupt:
    # except (KeyboardInterrupt, SystemExit):
        print(": ctrl-c event")
        observer.stop()
        observer.join()
        os._exit(0)
        # print("catching a ctrl-c event")
        # # my_client.stop()
        # die_nicely(oper)
