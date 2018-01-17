"""
File Routines for handling fastq files and monitoring locations. Built on watchdog.
"""
import os,sys
import threading
import time
import gzip

from tqdm import tqdm
from apifunc import Runcollection
from Bio import SeqIO
from watchdog.events import FileSystemEventHandler


def parsefastq(fastq, rundict, args, header):
    #print ('processing reads')
    if fastq.endswith(".gz"):
        #we have a potentiall gzipped file so:
        with gzip.open(fastq, "rt") as handle:
            for record in tqdm(SeqIO.parse(handle, "fastq")):
                descriptiondict = parsedescription(record.description)
                if descriptiondict["runid"] not in rundict:
                    rundict[descriptiondict["runid"]] = Runcollection(args,header)
                    rundict[descriptiondict["runid"]].add_run(descriptiondict)
                rundict[descriptiondict["runid"]].add_read(record, descriptiondict, fastq)
    else:
        for record in tqdm(SeqIO.parse(fastq, "fastq")):
            descriptiondict = parsedescription(record.description)
            if descriptiondict["runid"] not in rundict:
                rundict[descriptiondict["runid"]] = Runcollection(args, header)
                rundict[descriptiondict["runid"]].add_run(descriptiondict)
            rundict[descriptiondict["runid"]].add_read(record, descriptiondict,fastq)
    for runs in rundict:
        rundict[runs].commit_reads()



def parsedescription(description):
    descriptiondict = dict()
    descriptors = description.split(" ")
    del descriptors[0]
    for item in descriptors:
        bits = item.split("=")
        descriptiondict[bits[0]] = bits[1]
    return descriptiondict


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
                if (f.endswith(".fastq") or f.endswith(".fastq.gz")):
                    file_list_dict[os.path.join(path, f)] = os.stat(os.path.join(path, f)).st_mtime
                    # try:
                    #    file_descriptor = update_file_descriptor(os.path.join(path, f),file_descriptor)
                    # except:
                    #    pass
    print("processed %s files" % (counter))
    print("found %d existing fastq files to process first." % (len(file_list_dict)))
    return file_list_dict


class FastqHandler(FileSystemEventHandler):
    def __init__(self, args,header):
        """Collect information about files already in the folders"""
        self.file_descriptor = dict()
        self.args = args
        self.header=header
        # adding files to the file_descriptor is really slow - therefore lets skip that and only update the files when we want to basecall thread_number
        self.creates = file_dict_of_folder_simple(args.watchdir)
        self.processing = dict()
        self.running = True
        self.rundict = dict()
        t = threading.Thread(target=self.processfiles)

        try:
            t.start()
            print("Watchdog started")
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
            starttime = time.time()
            print("processfiles running {}".format(starttime))
            for fastqfile, createtime in tqdm(sorted(self.creates.items(), key=lambda x: x[1])):
                delaytime = 0
                if (int(
                        createtime) + delaytime < time.time()):  # file created 5 sec ago, so should be complete. For simulations we make the time longer.
                    #print(fastqfile)
                    del self.creates[fastqfile]
                    parsefastq(fastqfile, self.rundict,self.args, self.header)
            print ("processed in {}".format(time.time()-starttime))
            readsuploaded=0
            for runid in self.rundict:
                print("RunID", runid)
                print("Read Number:", self.rundict[runid].readcount, "Total Length:", self.rundict[runid].cumulength,
                      "Average Length", self.rundict[runid].cumulength / self.rundict[runid].readcount, "Chan Count",
                      len(self.rundict[runid].chandict))
                readsuploaded+=self.rundict[runid].readcount
                (mean, median, std, maxval, minval) = self.rundict[runid].mean_median_std_max_min()
                print("mean", mean, "median", median, "std", std, "max", maxval, "min", minval)
                # print self.rundict[runid].timeid
                #self.rundict[runid].parse1minwin()
                #os._exit(0)
            print ("Uploaded {} reads in total.".format(readsuploaded))

            time.sleep(5)

    def on_created(self, event):
        """Watchdog counts a new file in a folder it is watching as a new file"""
        """This will add a file which is added to the watchfolder to the creates and the info file."""
        if (event.src_path.endswith(".fastq") or event.src_path.endswith(".fastq.gz")):
            self.creates[event.src_path] = time.time()

    def on_modified(self, event):
        if (event.src_path.endswith(".fastq") or event.src_path.endswith(".fastq.gz")):
            self.creates[event.src_path] = time.time()

    def on_moved(self, event):
        """Watchdog considers a file which is moved within its domain to be a move"""
        """When a file is moved, we just want to update its location in the master dictionary."""
        if (event.dest_path.endswith(".fastq") or event.dest_path.endswith(".fastq.gz")):
            print("seen a fastq file move")

    def on_deleted(self, event):
        print("On Deleted Called", event.src_path)
