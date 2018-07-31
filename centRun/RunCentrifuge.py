import subprocess
import os
from io import StringIO
import pandas as pd
from ete3 import NCBITaxa
import numpy as np
import threading
from queue import Queue
import time
import sys
import signal
from centRun.models import CentOutput, LineageKey, LineageValues, MetaGenomicsMeta
from jobs.models import JobMaster
from reads.models import FastqRead
from centRun.serializers import CentSerialiser
from collections import defaultdict


class Reader(threading.Thread):
    """Reader class for opening a thread, that constantly reads the lines that are appearing in the process.stdout
    from centrifuge, and adds it to a Queue. Inherits threading.Thread
    Params:
        out: The buffered reader for stdout from the subprocess pipe
        queue: The python queue that we are going to add to

    """

    def __init__(self, out, queue):
        threading.Thread.__init__(self)
        self.out = out
        self.queue = queue
        self.running = True

    def run(self):
        """
        When called, repeatedly flushes stdout and adds any lines to a queue
        :return:
        """
        while self.running:
            # sys.stdout.flush()
            self.out.flush()
            time.sleep(.5)
            for line in self.out:
                self.queue.put(line)

    def close(self):
        """
        When called shuts the while loop in Reader.run and closes the thread with sys.exit
        :return:
        """
        print("exiting thread")
        self.running = False
        sys.exit()


class Centrifuger(Reader):
    """
        The Centrifuger class is a collection of methods that are required to perform the open centrifuge analysis on
        batches of reads called from the database. Inherits the Reader class.
    """

    def __init__(self, flowcell_id, flowcell_job_id):
        """

        Parameters
        ----------
        flowcell_id: The id of the flowcell passed into the celery task. Organises centrifuge output objects in the database
        so they can be retrieved together, and so reads can be collected
        flowcell_job_id: The id the job is stored under in the jobmasters table
        """
        self.flowcell_id = flowcell_id
        # The id this job is stored under in the Jobmaster table
        self.flowcell_job_id = flowcell_job_id
        # make the queue object available to the whole class
        self.outputLinesQueue = None
        # Continue centrifuging
        self.scan = True
        # Thread available to the whole class
        self.thread1 = None
        # The number of the last read centrifuged
        self.last_read = 0
        self.limit = 0

    @staticmethod
    def delete_series(series, df):
        """

        Parameters
        ----------
        series : (type list of Series) - series to be deleted
        df : (type a pandas DataFrame) - DataFrame to delete them from

        :return : a DataFrame with the series in the list deleted
        -------

        """
        for s in series:
            if type(df[s]).__name__ == "Series":
                df = df.drop(s, 1)
        return df

    def start_thread(self, pipe):
        """

        Parameters
        ----------
        pipe: (type: bytestream pipe) - The pipe connected to the centrifuge process

        Returns - No thing
        -------

        """
        # Queue, a python object that is first in first out, and thread safe,
        # to store individual centrifuge output lines
        self.outputLinesQueue = Queue(maxsize=0)
        self.thread1 = Reader(pipe, self.outputLinesQueue)
        # dies with the program
        self.thread1.setDaemon(True)
        self.thread1.start()

    @staticmethod
    def start_centrifuge():
        print("centrifuge index loaded")
        cmd = "centrifuge -f --mm -x /home/rory/testData/p_compressed -"
        # capture stdOut and input is the created fastq encode in utf-8 in bytes
        process = subprocess.Popen(cmd, shell=True, preexec_fn=os.setsid)

    @property
    def run_centrifuge(self):
        """

        Returns nothing.
        -------
         Run the centrifuge command, and open the subprocess to run it. Keep the process open and pipe the reads into
         it Write the reads incrementally to the processes stdin, capture the stdout and parse it into a useful format.
         Once parsed it's deposited into the database, updating any exisiting records. The organising key is a uuid,
         stored in self.foreign key.
        """

        cmd = "centrifuge -f --mm -x /home/rory/testData/p_compressed -"
        # capture stdOut and input is the created fastq encode in utf-8 in bytes
        process = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1,
                                   preexec_fn=os.setsid)
        # Start the thread that captures the stdout and adds to the queue
        self.start_thread(process.stdout)
        # boolean, if you have the header line from centrifuge
        header_line_got = False
        # Initialise header string
        header = ""
        # Counter for totla centOut
        total_centout = 0
        # Instance of the Ncbi taxa class
        ncbi = NCBITaxa()
        # check the organising uuid
        print(self.flowcell_id)
        # Loop whilst self.scan = True
        while self.scan:
            print(f"id is {self.flowcell_id}")
            print(f"flowcell_job_id is {self.flowcell_job_id}")
            JobMaster.objects.filter(pk=self.flowcell_job_id).update(running=True)
            # return all currently present reads
            cursor = FastqRead.objects.filter(run__flowcell_id__in={self.flowcell_id})
            print(cursor.values())
            doc_no = len(cursor)
            self.limit = doc_no - self.last_read
            last_read = list(cursor)[-1].id
            # Get the metadata for this run
            metadata = MetaGenomicsMeta.objects.get(meta_id=self.foreign_key)
            # The number of reads that we have in the database
            metadata.number_of_reads = doc_no
            metadata.save()
            print(f"lastread is {last_read}")
            print(f"self.last read is {self.last_read}")
            # if last read is too close to new last read
            if last_read - self.last_read < 100:
                self.scan = False
                print("If")
            elif self.last_read == last_read:
                print("Elif")
                self.scan = False
                break
            else:
                print("Else")
                self.last_read = last_read
            print(f"found {doc_no} reads in the database")

            # create fastq string by joining the read_id and sequence in the format, for all docs in cursor
            fastq = "".join([str(">" + c.read_id + "\n" + c.fastqreadextra.sequence + "\n") for c in cursor])

            # run the centrifuge command, the ending dash is EXTREMELY important
            print("centrifuging...")

            # Write the generated fastq file to stdin, passing it to the command
            process.stdin.write(fastq.encode("utf-8"))
            process.stdin.flush()
            print("centrifuged")
            # get all items in the queue
            q_size = self.outputLinesQueue.qsize()
            # a list to store the elements (Centrifuge output lines) as they are removed fromm the queue
            queue_out = []
            # if you have already got the first line of the output
            if header_line_got:
                # make sure the header line from the centrifuge output is first at the list each time
                queue_out.insert(0, header)
            # iterate over the loop and get items from the queue
            for i in range(q_size):
                line = str(self.outputLinesQueue.get_nowait(), "utf-8")

                # if you don't have the header line, it's the very first line, get it and store it
                if not header_line_got:
                    header = line
                    header_line_got = True
                # append the line removed from the queue to an array
                queue_out.append(line)
                # tell the queue that that element has been dealt with
                self.outputLinesQueue.task_done()
            # create the cent out string (tab separated and pass to pandas read csv)
            # centOut is the stdOut output of the centrifuge process decoded into a stringIo object, is tab delimited
            cent_out = "".join(queue_out)
            # total number of lines of centrifuge output dealt with
            total_centout += q_size
            print(f"number of centrifuge output lines is {total_centout}")
            # output fields is the column headers for the pandas data frame
            output_fields = ["readID", "seqID", "taxID", "numMatches"]

            # create the DataFrame from the output
            df = pd.read_csv(StringIO(cent_out), sep="\t", usecols=output_fields)
            # create DataFrame column unique, where if the number of matches is 1, unique is 1, if matches is more than
            # 1 unique is 0, using numpy condition
            df["unique"] = np.where(
                df["numMatches"] == 1,  # condition
                1,  # True
                0  # False
            )
            # get all the taxIds from the dataframe
            taxid_list = df["taxID"].values

            # use the ncbi thing to get the species names for the tax ids in the tax ids list
            taxid2name = ncbi.get_taxid_translator(taxid_list)

            def get_taxa(tax_id):
                """
                if these lineages contain one of these taxIds set a new series in the dataframe row listing it's
                10239 -> Viruses
                2 -> bacteria
                4751 -> Fungi
                33208 -> Anamalia
                9605 -> Homo

                :param tax_id:
                :return: return the names of taxons that are in the lineage that we are interested in
                Parameters
                ----------
                taxID (type Int) the taxID to query for the lineage

                Returns - return taxa level string kingdom in
                -------

                """
                # taxId of 0 is unclassified
                if tax_id is 0:
                    return 0
                # get the lineage corresponding to that taxID
                taxid2lineage = ncbi.get_lineage(tax_id)
                # list of taxIds to look for in the lineage of the taxId in the dataframe row
                target_taxa = [10239, 2, 4751, 33208, 9605, 2157]
                # return the taxon, of taxIDS that have a lineage taxID that is in the list b
                taxons = list(filter(lambda x: x in target_taxa, taxid2lineage))
                # get the name of thhe taxiD, as it's adict get the values
                taxons = ncbi.get_taxid_translator(taxons).values()
                # if lineage value is unclassified
                if not taxons:
                    return "Alex"
                else:
                    return ",".join(str(s) for s in taxons)
            # map over that series taXID to create a new series for the taxa
            # df["taxa"] = df["taxID"].map(get_taxa)

            # insert the taxid2name dict items into the dataframe name columns
            df["name"] = df["taxID"].map(taxid2name)
            # any NaNs replace with Unclassified
            df["name"] = df["name"].fillna("Rory")
            #  set the index to be taxID
            df = df.set_index("taxID")
            # create lots of mini data frames grouping by taxID
            gb = df.groupby(level="taxID")
            # create a sumUnique column in the DataFrame, where the value
            # is the sum of the unique column in the corresponding
            # grouped_by table
            df["sumUnique"] = gb["unique"].sum()
            # create a number of reads column in the data frame by getting
            # the number of rows in the corresponding grouped_by
            # table
            df["numReads"] = gb.size()
            # delete these columns, no longer needed
            df = self.delete_series(["readID", "seqID", "numMatches", "unique"], df)

            # remove duplicate rows in the data frame,so we only have one entry for each species
            df = df.drop_duplicates(keep="first")
            # rename columns to match database models
            # print(df.head(n=10))
            # if the number of documents is the same as the limit keep the same skip for the next query
            # print(f"limit is {self.limit}")
            # print(f"skip was {self.skip}")
            # if theres more than the number of skip left above the limit, add another skip to the milit
            file_lines = fastq.count("\n")
            # else set the skip to however many documents you just pulled from the database, i.e the remainder

            self.last_read += file_lines / 2
            int(self.last_read)

            print(f"last read is now {self.last_read}")

            # get the previous entries in the database
            # =========================================== # database dataframe
            # reset index
            # query the index get all the objects in the database at this point
            queryset = CentOutput.objects.filter(task_meta=self.flowcell_id).values()
            print(f"the length of the queryset is {len(queryset)}")
            if queryset:
                print("there was a awueryset")
                df = df.reset_index()
                # create a dictionary from queryset with taxID as key
                prev_entries_dict = {int(x["taxID"]): x for x in queryset}
                # create dataframe from prev_entries_dict fill from database
                prev_df = pd.DataFrame.from_dict(prev_entries_dict, orient="index")
                # reset index on dataframe
                prev_df = prev_df.reset_index()
                # rename columns to match the dataframe columns in dataframe
                #  generated from the centrifuge query this iteration
                prev_df = prev_df.rename(columns={'num_reads': 'numReads', 'sum_unique': 'sumUnique'})
                # append the two dataframes together
                df = df.append(prev_df, ignore_index=True)
                # set taxID as index
                df = df.set_index("taxID")
                # group by taxID
                gb = df.groupby(level="taxID")
                # create an  updatedSumUnique column in the DataFrame, where the value
                # is the sum of the unique column in the corresponding
                # grouped_by table
                df["updatedSumUnique"] = gb["sumUnique"].sum()
                # create a number of updatedNumReads column in the data frame by getting
                # the number of rows in the corresponding grouped_by table, for total in each group
                df["updatedNumReads"] = gb["numReads"].sum()
                print("\n")
                print(df.head(15))
                # delete the unnecessary old values series
                df = self.delete_series(["numReads", "sumUnique"], df)
                # rename updated to columns to names that model is expecting
                df = df.rename(columns={'updatedSumUnique': 'sum_unique', 'updatedNumReads': 'num_reads'})
                # reset index
                df = df.reset_index()
                # delete index series
                df = self.delete_series(["index"], df)
            else:
                df = df.reset_index()
                df = df.rename(columns={'sumUnique': 'sum_unique', 'numReads': 'num_reads'})
                # drop all duplicated lines to keep only one for each entry, Atomicity
            df = df.drop_duplicates(keep="first")
            # create bulk list, to populated with CentOutput items to insert
            bulk_list = []
            # create new defaultDict which creates a default dict for missing values, don't think we ever need this
            # default behaviour TODO double check this
            new = defaultdict(lambda: defaultdict())
            # Get the taxIDs in the Dataframe
            taxid_list2 = df["taxID"].values
            # Get the taxIDs of lineages that are in the database, so we can use those
            already_lineaged_tax_ids = LineageValues.objects.values_list("tax_id", flat=True)
            # Check the taxIDs we have in the dataframe, to see if they have matching lineages already.
            # If not the set subtraction means that unique, lineageless taxIDs end up in the not_prev_lineaged dict
            not_prev_lineaged = list(set(taxid_list2) - set(already_lineaged_tax_ids))
            # Get the lineages for these taxIDs in the list, returns a dict keyed by taxID,
            # value is listof taxIDS in lineage
            lineages_taxidlist_dict = ncbi.get_lineage_translator(not_prev_lineaged)
            # filter out no ranks, like root and cellular organism
            tax_rank_filter = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

            print("Entering lineages loop")
            # loop over the new lineages dictionary, the values of which is a list of taxIDs from root to the organism,
            # key is a int taxID
            for key, value in lineages_taxidlist_dict.items():
                # get the ranks of a list of taxID. Returns a dict with keys of the taxID in the list of taxIDS,
                # value is the rank, i.e root, kingdom etc.
                lineage_ranked = ncbi.get_rank(value)
                print("taxid_species_lookup_dict")
                # a dict where key is taxID, value is the name i.e {562: Escherichia coli}
                taxid_species_lookup_dict = ncbi.get_taxid_translator(value)
                print("new")
                # Populate the defaultDict, by looping the lineage ranked dict and adding it where they fit the
                # condition, the if key==taxID is if they are subspecies or strain
                new[key] = {rank: taxid_species_lookup_dict[tax_id] for tax_id, rank in lineage_ranked.items()
                            if rank in tax_rank_filter or key == tax_id}
            print("lindf")
            # create dataframe, the .T means that the dataframe is transposed, so dict key is the series key
            lin_df = pd.DataFrame.from_records(new).T
            # rename the no rank column to strain
            lin_df = lin_df.rename(columns={"no rank": "strain"})
            # if these series don't exist add them and populate with numpy NaN so there's no key errors
            for series in ["subspecies", "strain", "subStrainSpecies"]:
                if series not in lin_df.keys():
                    lin_df[series] = np.NaN
            # some strains are strains of a subspecies. But ncbi only returns the terminal (leaf) which is the
            # strain and not the subspecies above it so this function splits out the subspecies and adds it to a new
            # series

            def subspecies_determine(name):
                """

                Parameters
                ----------
                name - (type string) The name of the strain

                Returns - (type String) The subspecies as a string
                -------

                """
                # if the name is a string, not a float
                if type(name) is str:

                    if "subsp." in name:
                        print(name)
                        basestring = name.split("subsp.", 1)[0]
                        string = name.split("subsp.", 1)[1]
                        string = string.split(" ")[1]
                        subspecies = basestring + " subsp. " + string
                        print(subspecies)
                        return subspecies
                else:
                    pass
            # lin_df.fillna("NA")
            print("determining subspecies")
            # create new additional subspecies column, from strains column which has subspecies above it in taxa level
            lin_df["subStrainSpecies"] = lin_df["strain"].map(subspecies_determine)
            # merge new subspecies column with existing column
            lin_df["subspecies"].fillna(lin_df["subStrainSpecies"], inplace=True)
            # delete the new subspecies column
            self.delete_series(["subStrainSpecies"], lin_df)
            insert_list = []
            print("entering for first append")
            # iterate over lineage dataframe to save to database #TODO maybe rewrite into using mysql
            for ind, row in lin_df.iterrows():
                print(int(ind))
                # Add taxID to linage key model, where already lineaged taxIDs are store for fast lookup
                fk = LineageKey(tax_id=int(ind))
                fk.save()
                print("saving")
                # add to insert list the created lineage value objects to bulk create them
                insert_list.append(LineageValues(superkingdom=row["superkingdom"], phylum=row["phylum"],
                                                 tax_id=int(ind), classy=row["class"], order=row["order"],
                                                 family=row["family"], genus=row["genus"], species=row["species"],
                                                 subspecies=row["subspecies"], strain=row["strain"]
                                                 ))
            print("bulk creating")
            LineageValues.objects.bulk_create(insert_list)
            # entering second thing
            print("entering second thing")
            print(df.head(n=15))
            # iterate over the centrifuge output dataframe
            for ind, row in df.iterrows():
                # create CentOutput object #TODO if else, that will allow skipping of serialisation below
                insert = CentOutput(name=row["name"], taxID=row["taxID"],
                                    num_reads=row["num_reads"], sum_unique=row["sum_unique"],
                                    task_meta=self.flowcell_id, flowcell_id=self.flowcell_id)

                bulk_list.append(insert)
            # If there is already objects in teh database, do this update or create
            if CentOutput.objects.filter(task_meta=self.flowcell_id):
                # iterate create or update for everything that we are inserting
                for b in bulk_list:
                    # get the data from the queryset object
                    b_in = CentSerialiser(b)
                    CentOutput.objects.update_or_create(
                        taxID=b_in.data["taxID"], task_meta=self.flowcell_id, defaults=b_in.data
                    )

            # Bulk insert new model objects using bulk_create
            else:
                CentOutput.objects.bulk_create(bulk_list)
                print("Inserted first time")
            # if all reads are centrifuged
            JobMaster.objects.filter(pk=self.flowcell_job_id).update(last_read=last_read)
            metadata.reads_classified = doc_no
            metadata.save()
            # TODO need way to stop when all reads are done.
            self.scan = False
            if not self.scan:
                # stop the scan while loop
                JobMaster.objects.filter(pk=self.flowcell_job_id).update(running=False, complete=True)
                # set running to false. This means client stops querying
                # MetaGenomicsMeta.objects.filter(meta_id=self.foreign_key).update(running=False)
                print("Finished all reads in database")

                # kill centrifuge
                os.killpg(os.getpgid(process.pid), signal.SIGTERM)
                # close the thread
                self.thread1.close()
                self.thread1.join(timeout=1)
        print("finished")
        print(f"total centOut lines {total_centout}")
