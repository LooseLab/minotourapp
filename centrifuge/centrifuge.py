import subprocess
from io import StringIO
import pandas as pd
from ete3 import NCBITaxa
import numpy as np
from centrifuge.models import CentOutput, LineageValues, MetaGenomicsMeta, SankeyLinks
from jobs.models import JobMaster
from reads.models import FastqRead
from collections import defaultdict
from celery import task
from django.utils import timezone
from datetime import datetime
import time as timmy
from minotourapp.utils import get_env_variable


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
        # Continue centrifuging
        self.scan = True
        # The number of the last read centrifuged
        self.last_read = 0
        # The number of reads to skip when getting the reads out of the database
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
        centrifuge_path = get_env_variable("MT_CENTRIFUGE")
        index_path = get_env_variable("MT_CENT_INDEX")
        cmd = "perl" + centrifuge_path + " -f --mm -x " + index_path + " -"
        # Counter for total centOut
        total_centout = 0
        # Instance of the Ncbi taxa class, for taxonomic id manipulation
        ncbi = NCBITaxa()
        # Loop whilst self.scan = True
        d = datetime.now()
        # date is great
        # time is lime
        time = "{:%H:%M:%S}".format(d)
        # Get the task record from the JobMaster table in the database. Used to separate the results inserted,
        # and tell the Javascript and RunMonitor the task is running and when it is complete
        job_master = JobMaster.objects.get(pk=self.flowcell_job_id)
        job_master.running = True
        job_master.save()
        # Create a MetaData object about the Classification analysis, used to populate
        # the header on the visualisation.html page
        MetaGenomicsMeta(run_time=time, flowcell_id=self.flowcell_id, running=True, number_of_reads=0,
                         reads_classified=0, task=job_master).save()
        # Get the object that we just created back out
        metadata = MetaGenomicsMeta.objects.get(flowcell_id=self.flowcell_id, task__id=job_master.id)

        # While self.scan is true we query the fastqreads model for new readss that have appearedsince last time
        iteration_count = 0
        while self.scan:
            iteration_count += 1
            print("id is {}".format(self.flowcell_id))
            cursor = FastqRead.objects.filter(run__flowcell_id__in={self.flowcell_id})
            # return all currently present reads
            # Get the total number of reads in the database for this flowcell
            if not cursor:
                print("no reads in database")
                timmy.sleep(10)
                continue
            doc_no = int(cursor.count())
            # Update the read count on the JobMaster
            job_master.read_count = doc_no
            job_master.save()
            # Get the Id of the last read that we retrieve
            last_read_id = doc_no

            # The number of reads that we have in the database
            metadata.number_of_reads = doc_no
            metadata.save()

            print("found {} reads in the database for metagenomics on flowcell {}".format(doc_no, self.flowcell_id))
            print("last_read_id is {} and self.last_read_id is {}".format(last_read_id, self.last_read))

            # If the last read id that we have retrieved is the same as it was last iterations
            if self.last_read == last_read_id:
                # set the finish time
                print("if all reads finished")
                self.scan = False
                break
            # If less than 500 reads have appeared since last time, centrifuge them then stop scanning
            elif last_read_id - self.last_read < 100:
                print("elif less than 100")
                self.scan = False
            # Else Just update the last reads id on the self object and continue
            else:
                print("elsey more reads")
                self.last_read = last_read_id

            # Get all the reads skipping all we did last time

            print("slef.skip is {} and limit is {}".format(self.skip, doc_no))
            # create python list for zipping
            sequence_data = list(cursor.values_list("fastqreadextra__sequence", flat=True)[self.skip:doc_no])
            # create list of read ids for zipping
            read_ids = list(cursor.values_list("read_id", flat=True)[self.skip:doc_no])
            # Create list of tuples where
            # the 1st element is the read_id and the second is the sequence
            tupley_list = list(zip(read_ids, sequence_data))
            # update skip number
            self.skip = doc_no
            print("\n \n \n")
            print("number of reads in this iteration is {}".format(len(tupley_list)))
            print("\n \n \n")
            print("\n \n \n")
            print(tupley_list[-1])
            print("\n \n \n")
            # create fastq string by joining the read_id and sequence in the format, for all docs in cursor
            # fastq_list = []
            # for c in tupley_list:
            #     print(c)
            #     if c[1] is not None:
            #         fastq_list.append(str(">" + c[0] + "\n" + c[1] + "\n"))
            # fastq = "".join(fastq_list)
            fastq = "".join([str(">" + c[0] + "\n" + c[1] + "\n") for c in tupley_list if c[1] is not None])
            # Remove large objects to free up memory
            del sequence_data
            del read_ids
            del tupley_list
            # Write the generated fastq file to stdin, passing it to the command
            print("centrifuging...")
            # Use Popen to run the centrifuge command
            out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE) \
                .communicate(input=str.encode(fastq))
            # Remove large memory occupying string
            del fastq
            # the number of centrifuge lines
            cent_out = out.decode()
            # total number of lines of centrifuge output dealt with
            total_centout += cent_out.count("\n") - 1
            print("number of centrifuge output lines is {}".format(total_centout))
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
            # get np array of the taxids
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

            # =========================================== database dataframe PREVIOUS RESULTS
            # query the index get all the objects in the database at this point

            queryset = CentOutput.objects.filter(flowcell_id=self.flowcell_id, task__id=job_master.id).values()

            # Num matches series to map onto sankey links df below
            new_cent_num_matches = df["num_matches"]
            # Get the tax ids from the dataframe containg the newly produced centrifuge results
            new_tax_ids_list = df.index.values
            # filter out no ranks, like root and cellular organism
            tax_rank_filter = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
            # if theres a new dataset
            if queryset:
                df = df.reset_index()
                prev_df = pd.DataFrame(list(queryset))
                prev_df = self.delete_series(["id", "flowcell_id", "task_id"], prev_df)

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
                # delete the unnecessary old values series
                df = self.delete_series(["num_matches", "sum_unique"], df)
                # rename updated to columns to names that model is expecting
                df = df.rename(columns={'updated_sum_unique': 'sum_unique', 'updated_num_matches': 'num_matches'})
                # reset index
                df = df.reset_index()
                # delete index series
                # df = self.delete_series(["index"], df)
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
            # If not the set subtraction means that unique, lineage-less taxIDs end up in the not_prev_lineaged dict
            not_prev_lineaged = list(set(taxid_list2) - set(already_lineaged_tax_ids))

            # Get the lineages for these taxIDs in the new centrifuge output,
            # returns a dict keyed by taxID,
            # value is listof taxIDS in lineage
            lineages_taxidlist_dict = ncbi.get_lineage_translator(new_tax_ids_list)

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
            lin_df.head(n=20)
            # Order the DataFrame descending hierarchically, from kingdom to species
            sankey_lin_df = lin_df[tax_rank_filter]
            # Baqckfill the dataframe nas with the lowest common entry, if family is missing, and genus
            #  isn't, family NaN becomes the genus entry
            sankey_lin_df.fillna(axis=1, method="bfill", inplace=True)
            # Add the number of reads of the lowest common ancestor of a row to the dataframe
            # TODO This is the line that gives the warning on setting value on copy of slice
            sankey_lin_df["num_matches"] = new_cent_num_matches
            # initialise the series that will become our sankey pandas DataFrame
            source = pd.Series()
            target = pd.Series()
            value = pd.Series()
            tax_id = pd.Series()
            # rank = pd.Series()
            target_tax_level = pd.Series()
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
                    # rank = rank.append(sankey_lin_df["rank"])
                    sankey_lin_df["target_tax_level"] = target_clade
                    target_tax_level = target_tax_level.append(sankey_lin_df["target_tax_level"])
            # Create a DataFrame of links
            source_target_df = pd.concat([source, target, value, tax_id, target_tax_level], axis=1)
            # Rename the columns
            source_target_df.columns = ["source", "target", "value", "tax_id", "target_tax_level"]
            # Create a flowcell id series so each link has a flowcell id value
            source_target_df["flowcell_id"] = self.flowcell_id
            # Taxids from the current results that are in this df
            current_links_taxids = set(source_target_df["tax_id"])
            # TAxids that have links from this analysis already
            prev_links_taxids = set(SankeyLinks.objects.filter(flowcell_id=self.flowcell_id, task__id=job_master.id)
                                    .values_list("tax_id", flat=True))
            # Not yet linked is the taxIDs that don't have any result in the database
            not_yet_linked = current_links_taxids - prev_links_taxids
            # Create a subset dataframe that contains only links for new taxIDS
            to_create_df = source_target_df[source_target_df["tax_id"].isin(not_yet_linked)]
            # Create  a new series that contains the task record for this in each row
            to_create_df["job_master"] = job_master

            # Create a subset df containing links that already have reuslts
            to_update_df = source_target_df[source_target_df["tax_id"].isin(prev_links_taxids)]
            to_update_df["job_master"] = job_master
            # Link to contain objects for bulk insertion
            sankey_link_insert_list = []
            if not to_create_df.empty:
                def sankey_bulk_insert_list(row):
                    """
                    apply to the dataframe and create sankeylinks objects for each row
                    :param row: The dataframe row
                    :return:
                    """
                    sankey_link_insert_list.append(SankeyLinks(flowcell_id=row["flowcell_id"], source=row["source"],
                                                               target=row["target"], value=row["value"],
                                                               tax_id=row["tax_id"],
                                                               task=row["job_master"],
                                                               target_tax_level=row["target_tax_level"]))
                    return
                # Perform the apply
                to_create_df.apply(sankey_bulk_insert_list, axis=1)
                # Bulk create the objects
                SankeyLinks.objects.bulk_create(sankey_link_insert_list)

            # ## UPDATE existing links ## #
            if not to_update_df.empty:
                # Create a dataframe from the results that are already in the database
                prev_links_df = pd.DataFrame(list(SankeyLinks.objects.filter(flowcell_id=self.flowcell_id,
                                                                             task__id=job_master.id).values()))
                # Create a subset df that contains already linked taxids
                to_combine_values_df = prev_links_df[prev_links_df["tax_id"].isin(prev_links_taxids)]
                # delete the redundant dataframe
                del prev_links_df

                to_update_df.set_index(["tax_id", "target_tax_level"], inplace=True)
                to_combine_values_df.set_index(["tax_id", "target_tax_level"], inplace=True)

                # Update the value by combining the series from the two dataframes
                to_update_df["value"] = to_update_df["value"] + to_combine_values_df["value"]

                to_update_df.reset_index(inplace=True)
                # delete unnecessary dataframe
                del to_combine_values_df

                def sankey_update_links(row):
                    """
                    Fetches the object and updates the value and rank
                    :param row: The dataframe row
                    :return:
                    """
                    SankeyLinks.objects.filter(tax_id=row["tax_id"], flowcell_id=row["flowcell_id"], source=row["source"],
                                               target=row["target"], task=row["job_master"])\
                        .update(value=row["value"])
                to_update_df.apply(sankey_update_links, axis=1)

            # ##################### Back to the non sankey calculations ####################
            # Create a dataframe for lineages that aren't already in the table
            new_lineages_for_saving_df = lin_df[lin_df.index.isin(not_prev_lineaged)]
            # if the df isn't empty
            if not new_lineages_for_saving_df.empty:
                # iterate over lineage dataframe to save to database #TODO rewrite into using mysql
                # set the index name so upon reset the added series has a name
                new_lineages_for_saving_df.index.name = "tax_id"
                new_lineages_for_saving_df.reset_index(inplace=True)
                # list to contain the bulk insert object
                lineages_insert_list = []

                def lineages_bulk_insert_list(row):
                    """
                    Apply function to the dataframe to poulate a list with a LineageValue for each row
                    :param row: The row of the dataframe
                    :return:
                    """
                    lineages_insert_list.append(LineageValues(superkingdom=row["superkingdom"], phylum=row["phylum"],
                                                              tax_id=row["tax_id"], classy=row["class"],
                                                              order=row["order"],
                                                              family=row["family"], genus=row["genus"],
                                                              species=row["species"],
                                                              subspecies=row["subspecies"], strain=row["strain"]
                                                              ))
                    return

                new_lineages_for_saving_df.apply(lineages_bulk_insert_list, axis=1)
                # Bulk create
                LineageValues.objects.bulk_create(lineages_insert_list)

            # List to contain CentOutput objects
            centoutput_insert_list = []
            # Give each row a flowcell id cell, as part of the flowcell_id series
            df["flowcell_id"] = self.flowcell_id
            # same with task from Jobmaster
            df["task"] = job_master

            def centoutput_bulk_list(row):
                """
                Append a CentoutPut object to a list  for each row in a datframe
                :param row: the row from the dataframe
                :return:
                """
                centoutput_insert_list.append(CentOutput(name=row["name"], tax_id=row["tax_id"],
                                                         num_matches=row["num_matches"], sum_unique=row["sum_unique"],
                                                         flowcell_id=row["flowcell_id"], task=row["task"]))
                return
            # apply the function to the rows of this dataframe
            df.apply(centoutput_bulk_list, axis=1)

            # If there is already objects in teh database, do this update or create
            if CentOutput.objects.filter(flowcell_id=self.flowcell_id, task__id=job_master.id):
                # iterate create or update for everything that we are inserting#
                for centoutput in centoutput_insert_list:

                    CentOutput.objects.update_or_create(
                        tax_id=centoutput.tax_id, flowcell_id=self.flowcell_id, task__id=job_master.id,
                        defaults={"name": centoutput.name, "tax_id": centoutput.tax_id,
                                  "num_matches": centoutput.num_matches, "sum_unique": centoutput.sum_unique,
                                  "flowcell_id": centoutput.flowcell_id, "task": centoutput.task}
                    )
            # Bulk insert new model objects using bulk_create
            else:
                CentOutput.objects.bulk_create(centoutput_insert_list)
                print("Inserted first time")
            # Update the jobmaster object fields that are relevant
            job_master.last_read = last_read_id
            job_master.save()
            metadata.reads_classified = doc_no
            metadata.save()
            if iteration_count < 20:
                print(iteration_count)
                timmy.sleep(4)
        # if self.scan is false
        if not self.scan:
            # Update the job_master object to running
            job_master.running = False
            job_master.complete = True
            job_master.save()
            # set running to false. This means client stops querying
            # MetaGenomicsMeta.objects.filter(meta_id=self.foreign_key).update(running=False)
            print("Finished all reads in database")

            start_time = metadata.timestamp.replace(tzinfo=None)
            end_time = timezone.now().replace(tzinfo=None)
            metadata.finish_time = str(end_time - start_time)
            metadata.save()
            print("finished")
            print("total centOut lines {}".format(total_centout))
