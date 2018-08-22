import subprocess
from io import StringIO
import pandas as pd
from ete3 import NCBITaxa
import numpy as np
from centrifuge.models import CentOutput, LineageValues, MetaGenomicsMeta, SankeyLinks
from jobs.models import JobMaster
from reads.models import FastqRead
from centrifuge.serializers import CentSerialiser
from collections import defaultdict
from celery import task
from django.utils import timezone

# TODO del unused df and things


class Centrifuger:
    """
        The Centrifuger class is a collection of methods that are required to perform the open centrifuge analysis on
        batches of reads called from the database. Inherits the Reader class.
    """

    def __init__(self, flowcell_id, flowcell_job_id):
        """

        Parameters
        ----------
        flowcell_id: The id of the flowcell passed into the celery task. Organises
        centrifuge output objects in the database
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
        self.skip = 0

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

    @property
    @task()
    def run_centrifuge(self):
        """

        Returns nothing.
        -------
         Run the centrifuge command, and open the subprocess to run it. Keep the process open and pipe the reads into
         it Write the reads incrementally to the processes stdin, capture the stdout and parse it into a useful format.
         Once parsed it's deposited into the database, updating any existing records. The organising key is a uuid,
         stored in self.foreign key.
        """

        cmd = "centrifuge -f --mm -x /home/rory/testData/p_compressed -"
        # Counter for total centOut
        total_centout = 0
        # Instance of the Ncbi taxa class, for taxonomic id manipulation
        ncbi = NCBITaxa()
        # Loop whilst self.scan = True

        job_master = JobMaster.objects.get(pk=self.flowcell_job_id)
        job_master.running = True
        job_master.save()
        metadata = MetaGenomicsMeta.objects.get(flowcell_id=self.flowcell_id, task__id=job_master.id)
        while self.scan:
            print(f"id is {self.flowcell_id}")
            # return all currently present reads
            cursor = FastqRead.objects.filter(run__flowcell_id__in={self.flowcell_id})
            doc_no = cursor.count()
            last_read_id = cursor[doc_no-1].id
            # Get the metadata for this run

            # The number of reads that we have in the database
            metadata.number_of_reads = doc_no
            metadata.save()
            print(f"lastread is {last_read_id}")
            print(f"self.last read is {self.last_read}")
            # if last read is too close to new last read
            print(f"last_read_id is {last_read_id} and self.last_read_id is {self.last_read}")

            if self.last_read == last_read_id:
                # set the finish time
                self.scan = False
                break
            elif last_read_id - self.last_read < 500:
                self.scan = False
            else:
                self.last_read = last_read_id
            print(f"found {doc_no} reads in the database")
            # Get all the reads skipping all we did last time
            cursor = FastqRead.objects.filter(run__flowcell_id__in={self.flowcell_id})[self.skip:doc_no]
            # create fastq string by joining the read_id and sequence in the format, for all docs in cursor
            # TODO sends individual queries to fetch the data for each read. Get the data in memory in one go?
            fastq = "".join([str(">" + c.read_id + "\n" + c.fastqreadextra.sequence + "\n") for c in cursor])

            # Write the generated fastq file to stdin, passing it to the command
            print("centrifuging...")
            out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE) \
                .communicate(input=str.encode(fastq))
            del fastq
            # the number of centrifuge lines
            cent_out = out.decode()
            # total number of lines of centrifuge output dealt with
            total_centout += cent_out.count("\n") - 1
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
            df.rename(columns={"taxID": "tax_id"}, inplace=True)
            taxid_list = df["tax_id"].values

            # use the ncbi thing to get the species names for the tax ids in the tax ids list
            taxid2name = ncbi.get_taxid_translator(taxid_list)

            # insert the taxid2name dict items into the dataframe name columns
            df["name"] = df["tax_id"].map(taxid2name)
            # any NaNs replace with Unclassified
            df["name"] = df["name"].fillna("Rory")
            #  set the index to be taxID
            df = df.set_index("tax_id")
            # create lots of mini data frames grouping by taxID
            gb = df.groupby(level="tax_id")
            # create a sumUnique column in the DataFrame, where the value
            # is the sum of the unique column in the corresponding
            # grouped_by table
            df["sum_unique"] = gb["unique"].sum()
            # create a number of reads column in the data frame by getting
            # the number of rows in the corresponding grouped_by
            # table
            df["num_matches"] = gb.size()
            # delete these columns, no longer needed
            df = self.delete_series(["readID", "seqID", "numMatches", "unique"], df)

            # remove duplicate rows in the data frame,so we only have one entry for each species
            df = df.drop_duplicates(keep="first")
            # rename columns to match database models
            # print(df.head(n=10))
            # if the number of documents is the same as the limit keep the same skip for the next query
            # if theres more than the number of skip left above the limit, add another skip to the milit

            # get the previous entries in the database
            # =========================================== database dataframe
            # reset index
            # query the index get all the objects in the database at this point
            queryset = CentOutput.objects.filter(flowcell_id=self.flowcell_id, task__id=job_master.id).values()
            print(f"the length of the queryset is {len(queryset)}")

            new_cent_num_matches = df["num_matches"]
            # Get the tax ids from the dataframe containg the newly produced centrifuge results
            new_tax_ids_list = df.index.values
            # filter out no ranks, like root and cellular organism
            tax_rank_filter = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

            if queryset:
                df = df.reset_index()

                prev_df = pd.DataFrame(list(queryset))
                # append the two dataframes together
                df = df.append(prev_df, ignore_index=True)
                # set taxID as index
                df = df.set_index("tax_id")
                # group by taxID
                gb = df.groupby(level="tax_id")
                # create an  updatedSumUnique column in the DataFrame, where the value
                # is the sum of the unique column in the corresponding
                # grouped_by table
                df["updated_sum_unique"] = gb["sum_unique"].sum()
                # create a number of updatednum_matches column in the data frame by getting
                # the number of rows in the corresponding grouped_by table, for total in each group
                df["updated_num_matches"] = gb["num_matches"].sum()
                print("\n")
                print(df.head(15))
                # delete the unnecessary old values series
                df = self.delete_series(["num_matches", "sumUnique"], df)
                # rename updated to columns to names that model is expecting
                df = df.rename(columns={'updated_sum_unique': 'sum_unique', 'updated_num_matches': 'num_matches'})
                # reset index
                df = df.reset_index()
                # delete index series
                df = self.delete_series(["index"], df)
            else:
                df = df.reset_index()
                # drop all duplicated lines to keep only one for each entry, Atomicity
            df = df.drop_duplicates(keep="first")

            # create new defaultDict which creates a default dict for missing values, don't think we ever need this
            # default behaviour
            new = defaultdict(lambda: defaultdict())
            # Get the taxIDs in the Dataframe
            taxid_list2 = df["tax_id"].values
            # Get the taxIDs of lineages that are in the database, so we can use those
            already_lineaged_tax_ids = LineageValues.objects.values_list("tax_id", flat=True)
            # Check the taxIDs we have in the dataframe, to see if they have matching lineages already.
            # If not the set subtraction means that u   nique, lineage-less taxIDs end up in the not_prev_lineaged dict
            not_prev_lineaged = list(set(taxid_list2) - set(already_lineaged_tax_ids))

            # Get the lineages for these taxIDs in the new centrifuge output,
            # returns a dict keyed by taxID,
            # value is listof taxIDS in lineage
            lineages_taxidlist_dict = ncbi.get_lineage_translator(new_tax_ids_list)

            print("Entering lineages loop")
            # loop over the new lineages dictionary, the values of which is a list of taxIDs from root to the organism,
            # key is a int taxID
            # TODO rewrite into pandas mapping function, low priority
            for key, value in lineages_taxidlist_dict.items():
                # get the ranks of a list of taxID. Returns a dict with keys of the taxID in the list of taxIDS,
                # value is the rank, i.e root, kingdom etc.
                lineage_ranked = ncbi.get_rank(value)
                # a dict where key is taxID, value is the name i.e {562: Escherichia coli}
                taxid_species_lookup_dict = ncbi.get_taxid_translator(value)
                # Populate the defaultDict, by looping the lineage ranked dict and adding it where they fit the
                # condition, the if key==taxID is if they are subspecies or strain
                new[key] = {rank: taxid_species_lookup_dict[tax_id] for tax_id, rank in lineage_ranked.items()
                            if rank in tax_rank_filter or key == tax_id}
            print("lindf")
            # create dataframe, the .T means that the dataframe is transposed, so dict key is the series title
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
            # ##############  Storing snakey links and nodes ################

            # Order the DataFrame descending hierarchically, from kingdom to species
            sankey_lin_df = lin_df[tax_rank_filter]
            # Baqckfill the dataframe nas with the lowest common entry, if family is missing, and genus
            #  isn't, family NaN becomes the genus entry
            sankey_lin_df.fillna(axis=1, method="bfill", inplace=True)
            # Add the number of reads of the lowest common ancestor of a row to the dataframe
            # TODO This is the line that gives the warning on setting value on copy of slice
            sankey_lin_df["num_matches"] = new_cent_num_matches
            # Sort the dataframe by the nummatches, descending
            sankey_lin_df.sort_values(["num_matches"], ascending=False, inplace=True)
            # Add a rank for the species
            sankey_lin_df["rank"] = np.arange(len(sankey_lin_df))
            # initialise the series that will become our sankey pandas DataFrame
            source = pd.Series()
            target = pd.Series()
            value = pd.Series()
            tax_id = pd.Series()
            rank = pd.Series()
            # Iterate across the dataframe columns, to create the source, target, num reads series
            for index, clade in enumerate(tax_rank_filter):
                source_clade = clade
                k = index + 1
                if k < len(tax_rank_filter):
                    target_clade = tax_rank_filter[k]
                    source = source.append(sankey_lin_df[source_clade])
                    target = target.append(sankey_lin_df[target_clade])
                    value = value.append(sankey_lin_df["num_matches"])
                    tax_id = tax_id.append(pd.Series(sankey_lin_df.index.values, index=sankey_lin_df.index))
                    rank = rank.append(sankey_lin_df["rank"])
            # Create a DataFrame of links
            source_target_df = pd.concat([source, target, value, tax_id, rank], axis=1)
            # Rename the columns
            source_target_df.columns = ["source", "target", "value", "tax_id", "rank"]
            # Set a multi index of source to target

            source_target_df["flowcell_id"] = self.flowcell_id

            current_links_taxids = set(source_target_df["tax_id"])

            prev_links_taxids = set(SankeyLinks.objects.filter(flowcell_id=self.flowcell_id, task__id=job_master.id)
                                    .values_list("tax_id", flat=True))

            not_yet_linked = current_links_taxids - prev_links_taxids

            to_create_df = source_target_df[source_target_df["tax_id"].isin(not_yet_linked)]

            to_create_df["job_master"] = job_master

            to_update_df = source_target_df[source_target_df["tax_id"].isin(prev_links_taxids)]

            sankey_link_insert_list = []

            if not to_create_df.empty:
                def sankey_bulk_insert_list(row):
                    sankey_link_insert_list.append(SankeyLinks(flowcell_id=row["flowcell_id"], source=row["source"],
                                                               target=row["target"], value=row["value"],
                                                               tax_id=row["tax_id"], rank=row["rank"],
                                                               task=row["job_master"]))
                    return

                to_create_df.apply(sankey_bulk_insert_list, axis=1)

                SankeyLinks.objects.bulk_create(sankey_link_insert_list)

            # ## UPDATE existing links ## #
            if not to_update_df.empty:
                prev_links_df = pd.DataFrame(list(SankeyLinks.objects.filter(flowcell_id=self.flowcell_id,
                                                                             task__id=job_master.id).values()))

                to_combine_values_df = prev_links_df[prev_links_df["tax_id"].isin(prev_links_taxids)]

                del prev_links_df

                to_combine_values_df.set_index("tax_id", inplace=True)

                to_update_df.set_index("tax_id", inplace=True)

                to_update_df["value"] = to_update_df["value"] + to_combine_values_df["value"]

                to_update_df.reset_index(inplace=True)

                del to_combine_values_df

                def sankey_update_links(row):
                    SankeyLinks.objects.filter(flowcell_id=row["flowcell_id"], source=row["source"],
                                               target=row["target"]).update(value=row["value"], rank=row["rank"])
                to_update_df.apply(sankey_update_links, axis=1)

            # ##################### Back to the non sankey calculations ####################

            new_lineages_for_saving_df = lin_df[lin_df.index.isin(not_prev_lineaged)]
            if not new_lineages_for_saving_df.empty:
                # iterate over lineage dataframe to save to database #TODO rewrite into using mysql
                new_lineages_for_saving_df.index.name = "tax_id"
                new_lineages_for_saving_df.reset_index(inplace=True)
                lineages_insert_list = []

                def lineages_bulk_insert_list(row):
                    lineages_insert_list.append(LineageValues(superkingdom=row["superkingdom"], phylum=row["phylum"],
                                                              tax_id=row["tax_id"], classy=row["class"],
                                                              order=row["order"],
                                                              family=row["family"], genus=row["genus"],
                                                              species=row["species"],
                                                              subspecies=row["subspecies"], strain=row["strain"]
                                                              ))
                    return

                new_lineages_for_saving_df.apply(lineages_bulk_insert_list, axis=1)
                LineageValues.objects.bulk_create(lineages_insert_list)

            # TODO rewrite iterrows to apply returnin a list
            centoutput_insert_list = []
            df["flowcell_id"] = self.flowcell_id
            df["task"] = job_master

            def centoutput_bulk_list(row):
                centoutput_insert_list.append(CentOutput(name=row["name"], tax_id=row["tax_id"],
                                                         num_matches=row["num_matches"], sum_unique=row["sum_unique"],
                                                         flowcell_id=row["flowcell_id"], task=row["task"]))
                return
            df.apply(centoutput_bulk_list, axis=1)

            # If there is already objects in teh database, do this update or create
            if CentOutput.objects.filter(flowcell_id=self.flowcell_id):
                # iterate create or update for everything that we are inserting
                for centoutput in centoutput_insert_list:
                    # get the data from the queryset object
                    cent_in = CentSerialiser(centoutput)
                    CentOutput.objects.update_or_create(
                        tax_id=cent_in.data["taxID"], flowcell_id=self.flowcell_id, defaults=cent_in.data
                    )

            # Bulk insert new model objects using bulk_create
            else:
                CentOutput.objects.bulk_create(centoutput_insert_list)
                print("Inserted first time")
            # if all reads are centrifuged
            job_master.last_read = last_read_id
            job_master.save()
            metadata.reads_classified = doc_no
            metadata.save()
            # TODO need way to stop when all reads are done.
        if not self.scan:
            # stop the scan while loop
            job_master.running = False
            job_master.complete = True
            job_master.save()
            # set running to false. This means client stops querying
            # MetaGenomicsMeta.objects.filter(meta_id=self.foreign_key).update(running=False)
            print("Finished all reads in database")

            start_time = metadata.timestamp.replace(tzinfo=None)
            end_time = timezone.now().replace(tzinfo=None)
            metadata.finish_time = end_time - start_time
            metadata.save()
            print("finished")
            print(f"total centOut lines {total_centout}")
