import subprocess
from io import StringIO
import pandas as pd
from celery.utils.log import get_task_logger
from ete3 import NCBITaxa
import numpy as np
from centrifuge.models import CentOutput, LineageValues, MetaGenomicsMeta, SankeyLinks, CentOutputBarcoded
from jobs.models import JobMaster
from reads.models import FastqRead
from reads.models import Flowcell
from collections import defaultdict
from celery import task
from django.utils import timezone
from datetime import datetime
import time as timmy
from minotourapp.utils import get_env_variable
from django.db.models import ObjectDoesNotExist
from centrifuge.mapping import Metamap

logger = get_task_logger(__name__)


# TODO del unused df and things

def barcode_calculations(gb_df, ho_df):
    """
    Return the main cent output data frame with the barcode only result concatenated on
    :param gb_df: barcode data frame group data frame, one for each barcode
    :param ho_df: the main results data frame
    :return: A new df with the result for barcode appended
    """

    gb_df_gb = gb_df.groupby("tax_id")
    # calculate total unique matches for this species in this barcode
    sum_unique = gb_df_gb["unique"].sum()
    # calculate the total matches fot this species in this barcode
    num_matches = gb_df_gb.size()
    gb_df.drop(columns=["readID", "seqID", "numMatches", "unique"], inplace=True)
    gb_df.drop_duplicates(subset=["barcode", "name"], keep="first", inplace=True)

    # combine the dataframes
    gb_ap = pd.concat([sum_unique, num_matches, gb_df["name"]], axis=1)
    ho_df = ho_df.append(gb_ap)
    return ho_df


def centoutput_bulk_list(row):
    """
    Append a CentOutput object to a list  for each row in a dataframe
    :param row: the row from the dataframe
    :return: The list of newly created objects
    """
    return CentOutput(name=row["name"], tax_id=row["tax_id"],
                      flowcell=row["flowcell"], task=row["task"])


def bulk_create_list(row, job_master, flowcell_id):
    """
    Create barcoded lists
    :param job_master: The job_master object for this task run
    :param row: The dataframe row
    :param flowcell_id: the id of the flowcell
    :return: The  newly created objects
    """
    try:
        cent = CentOutput.objects.get(tax_id=row["tax_id"], task__id=job_master.id, flowcell__id=flowcell_id)
        return CentOutputBarcoded(num_matches=row["num_matches"],
                                  sum_unique=row["sum_unique"],
                                  output=cent,
                                  tax_id=row["tax_id"],
                                  barcode=row["barcode"])

    except ObjectDoesNotExist:
        logger.info("<<<<<")
        logger.info("Cent Object doesn't exist")
        logger.info("No matching entry for tax_id {} and task_id {}".format(row["tax_id"], job_master.id))
        logger.info("<<<<<")


def update_bar_values(row, flowcell_job_id, flowcell_id):
    """

    :param flowcell_id: The flowcell_id for this task
    :param flowcell_job_id: The id for the job_master for this task run
    :param row: The data frame row
    :return: The list of newly created objects
    """
    CentOutputBarcoded.objects.filter(tax_id=row["tax_id"], barcode=row["barcode"],
                                      output__task__id=flowcell_job_id,
                                      output__flowcell__id=flowcell_id).update(
    num_matches=row["updated_num_matches"], sum_unique=row["updated_sum_unique"])


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


def sankey_bulk_insert_list(row):
    """
    apply to the dataframe and create sankeylinks objects for each row
    :param row: The data frame row
    :return: The list of created objects
    """
    return SankeyLinks(flowcell=row["flowcell"], source=row["source"],
                       target=row["target"],
                       tax_id=row["tax_id"],
                       task=row["job_master"],
                       target_tax_level=row["target_tax_level"], value=row["value"],
                       barcode=row["barcode"])


def sankey_bulk_bar_insert(row, flowcell_id, flowcell_job_id):
    """
    Create or update the barcode objects in the database
    :param row: The dataframe row
    :param flowcell_id: The primary key for the flowcell the task is being run on
    :param flowcell_job_id: The Pk of the task ID
    :return: Nothing
    """

    # try:
    SankeyLinks.objects.filter(flowcell__id=flowcell_id, task__id=flowcell_job_id,
                               target_tax_level=row["target_tax_level"],
                               barcode=row["barcode"], tax_id=row["tax_id"]).update(value=row["updated_value"])
    # except KeyError as e:
    #     logger.warning(e.__traceback__)
    #     logger.info("<<<<<")
    #     logger.info("SankeyLinks doesn't exist")
    #     logger.info("No matching entry for tax_id {} and target level {} in barcode {}".format(row["tax_id"],
    #                                                                                            row["target_tax_level"],
    #                                                                                            row["barcode"]))
    #     logger.info("<<<<<")


def lineages_bulk_insert_list(row):
    """
    Apply function to the dataframe to populate a list with a LineageValue for each row
    :param row: The row of the dataframe
    :return: The list of newly created objects
    """
    return LineageValues(superkingdom=row["superkingdom"], phylum=row["phylum"],
                         tax_id=row["tax_id"], classy=row["class"],
                         order=row["order"],
                         family=row["family"], genus=row["genus"],
                         species=row["species"],
                         subspecies=row["subspecies"], strain=row["strain"]
                         )


class Centrifuger:
    """
        The Centrifuger class is a collection of methods that are required to perform the open centrifuge analysis on
        batches of reads called from the database. Inherits the Reader class.
    """

    def __init__(self, flowcell_job_id):
        """

        Parameters
        ----------
        flowcell_id: The id of the flowcell passed into the celery task. Organises
        centrifuge output objects in the database
        so they can be retrieved together, and so reads can be collected
        flowcell_job_id: The id the job is stored under in the jobmasters table
        """
        job_master = JobMaster.objects.get(pk=flowcell_job_id)

        self.flowcell_id = job_master.flowcell.id
        # The id this job is stored under in the Jobmaster table
        self.flowcell_job_id = flowcell_job_id
        # Continue centrifuging
        self.scan = True  # Remove?
        # The number of the last read centrifuged
        self.last_read = job_master.last_read
        # The number of reads to skip when getting the reads out of the database
        # self.skip = 1000
        # the number to chunk in
        self.chunk = 500

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

    def run_centrifuge(self):
        """

        Returns nothing.
        -------
         Run the centrifuge command, and open the subprocess to run it. Keep the process open and pipe the reads into
         it Write the reads incrementally to the processes stdin, capture the stdout and parse it into a useful format.
         Once parsed it's deposited into the database, updating any existing records. The organising key is a uuid,
         stored in self.foreign key.
        """
        job_master = JobMaster.objects.get(pk=self.flowcell_job_id)
        job_master.running = True
        job_master.save()
        logger.info("Metagenomics Flowcell id: {} - Running Metagenomics on Flowcell {}, with JobMaster ID {}".format(
            self.flowcell_id, self.flowcell_id, job_master.id)
        )

        centrifuge_path = get_env_variable("MT_CENTRIFUGE")
        index_path = get_env_variable("MT_CENTRIFUGE_INDEX")
        cmd = "perl " + centrifuge_path + " -f --mm -x " + index_path + " -"

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

        # Create a MetaData object about the Classification analysis, used to populate
        # the header on the visualisation.html page
        flowcell = Flowcell.objects.get(pk=self.flowcell_id)

        metadata = MetaGenomicsMeta(
            run_time=time,
            flowcell=flowcell,
            running=True,
            number_of_reads=0,
            reads_classified=0,
            task=job_master
        )

        metadata.save()

        # Get the object that we just created back out
        # metadata = MetaGenomicsMeta.objects.get(flowcell__id=self.flowcell_id, task__id=job_master.id)

        # While self.scan is true we query the fastqreads model for new readss that have appearedsince last time
        iteration_count = 0
        cursor = FastqRead.objects.filter(run__flowcell_id__in={self.flowcell_id})

        iteration_count += 1
        logger.info("id is {}".format(self.flowcell_id))

        # return all currently present reads
        # Get the total number of reads in the database for this flowcell
        # if not cursor:
        #     logger.info("no reads in database")
        #     timmy.sleep(10)
        #     continue
        doc_no = int(cursor.count())
        # doc_no = 500
        # Update the read count on the JobMaster
        # job_master.read_count = doc_no
        # job_master.save()
        # Get the Id of the last read that we retrieve
        # if doc_no < (self.skip + self.chunk):
        #     limit = doc_no
        # else:
        #     limit = self.skip + self.chunk
        # last_read_id = limit

        # The number of reads that we have in the database
        metadata.number_of_reads = doc_no
        metadata.save()
        logger.info(" Metagenomics Flowcell id: {} - found {} reads in the database for metagenomics on flowcell {}".format(
            self.flowcell_id,
            doc_no,
            self.flowcell_id)
        )
        logger.info("Metagenomics Flowcell id: {} - self.last_read_id is {}".format(
            self.flowcell_id,
            self.last_read)
        )

        # If the last read id that we have retrieved is the same as it was last iterations
        # if self.last_read == last_read_id:
        #     # set the finish time
        #     logger.info("if all reads finished")
        #     self.scan = False
        #     break
        # # If less than 500 reads have appeared since last time, centrifuge them then stop scanning
        # elif last_read_id - self.last_read < 100:
        #     logger.info("elif less than 100")
        #     self.scan = False
        # # Else Just update the last reads id on the self object and continue
        # else:
        #     logger.info("elsey more reads")
        #     self.last_read = last_read_id

        # Get all the reads skipping all we did last time

        logger.info("Metagenomics Flowcell id: {} - Reads Found".format(
            self.flowcell_id,
            doc_no)
        )
        limit = self.last_read + self.chunk
        # create python list for zipping
        sequence_data = list(cursor.values_list("fastqreadextra__sequence", flat=True)[self.last_read:limit])
        # create list of read ids for zipping
        read_ids = list(cursor.values_list("read_id", flat=True)[self.last_read:limit])
        # Create list of tuples where
        barcodes = list(cursor.values_list("barcode_name", flat=True)[self.last_read:limit])
        # the 1st element is the read_id and the second is the sequence
        tupley_list = list(zip(read_ids, sequence_data, barcodes))
        # update skip number
        # TODO UNCOMMENT
        # self.skip = limit
        logger.info("\n \n \n")
        logger.info("Metagenomics Flowcell id: {} - number of reads in this iteration is {}".format(self.flowcell_id,
                                                                                                    len(tupley_list)))

        # create fastq string by joining the read_id and sequence in the format, for all docs in cursor
        fastq = "".join([str(">read_id=" + c[0] + ",barcode=" + c[2] + "\n" + c[1] + "\n")
                         for c in tupley_list if c[1] is not None])

        # Remove large objects to free up memory
        del sequence_data
        del read_ids
        del tupley_list
        # Write the generated fastq file to stdin, passing it to the command
        logger.info("Metagenomics Flowcell id: {} - Loading index and Centrifuging".format(self.flowcell_id,
                                                                                           ))
        # Use Popen to run the centrifuge command
        out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE) \
            .communicate(input=str.encode(fastq))
        # Remove large memory occupying string
        del fastq
        # the number of centrifuge lines
        cent_out = out.decode()
        # total number of lines of centrifuge output dealt with
        total_centout += cent_out.count("\n") - 1
        logger.info("Metagenomics Flowcell id: {} - number of centrifuge output lines is {}".format(self.flowcell_id,
                                                                                                    total_centout))
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
        df = pd.concat([df, df["readID"].str.split(",").str[1].str.split("=").str[1]], axis=1)
        df.columns = ["readID", "seqID", "tax_id", "numMatches", "unique", "barcode"]
        df = pd.concat([df, df["readID"].str.split(",").str[0].str.split("=").str[1]], axis=1)
        df.columns = ["readID", "seqID", "tax_id", "numMatches", "unique", "barcode", "read_id"]

        # get np array of the taxids
        taxid_list = df["tax_id"].values
        # use the ncbi thing to get the species names for the tax ids in the tax ids list
        taxid2name_df = pd.DataFrame.from_dict(ncbi.get_taxid_translator(taxid_list), orient="index",
                                               columns=["name"])
        # insert the taxid2name dict items into the dataframe name columns
        df = pd.merge(df, taxid2name_df, how="outer", left_on="tax_id", right_index=True)
        # any NaNs replace with Unclassified
        df["name"].fillna("Unclassified", inplace=True)

        #  set the index to be taxID
        df = df.set_index("tax_id")
        # The barcoding
        gb_bc = df.groupby(["barcode"])
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
        # TODO future mapping target
        temp_targets = ["Escherichia coli", "Bacillus cereus"]

        name_df = df[df["name"].isin(temp_targets)]
        # if not name_df.empty:
            # map the target reads
            # m = Metamap(name_df, self.flowcell_id, self.flowcell_job_id)
            # m.map_the_reads()
        barcode_df = pd.DataFrame()
        # TODO vectorise these bad boys
        logger.info("Metagenomics Flowcell id: {} - The dataframe shape is {}".format(self.flowcell_id, df.shape))
        unique_barcode = set(barcodes)
        if len(unique_barcode) > 1:
            logger.info("Metagenomics Flowcell id: {} - These are the barcodes in this set".format(self.flowcell_id))
            logger.info(unique_barcode)
            barcode_df = gb_bc.apply(barcode_calculations, barcode_df)
            barcode_df.reset_index(inplace=True)
            barcode_df.rename(columns={"unique": "sum_unique", 0: "num_matches"}, inplace=True)
            barcode_df.set_index("tax_id", inplace=True)
        # delete these columns, no longer needed
        df.drop(columns=["readID", "seqID", "numMatches", "unique", "barcode", "read_id"], inplace=True)
        # remove duplicate rows in the data frame,so we only have one entry for each species
        df.drop_duplicates(keep="first", inplace=True)

        # =========================================== database dataframe PREVIOUS RESULTS
        # Get the tax ids from the dataframe containg the newly produced centrifuge results
        new_tax_ids_list = df.index.values
        # query the index get all the objects in the database for this analysis at this point
        queryset = CentOutput.objects.filter(flowcell__id=self.flowcell_id, task__id=job_master.id)

        # Num matches series to map onto sankey links df below

        # filter out no ranks, like root and cellular organism
        tax_rank_filter = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

        # drop all duplicated lines to keep only one for each entry, Atomicity
        df.reset_index(inplace=True)
        df = df.drop_duplicates(keep="first")
        df["barcode"] = "All reads"
        prev_df_tax_ids = list(queryset.values_list("tax_id", flat=True))
        # ###### BULK CREATE CENTOUTPUT OBJECTS ########
        # Give each row a flowcell object cell, as part of the flowcell_id series
        df["flowcell"] = flowcell
        # same with task from Jobmaster
        df["task"] = job_master
        cent_to_create = df[~df["tax_id"].isin(prev_df_tax_ids)]
        # apply row wise
        centoutput_insert_list = cent_to_create.apply(centoutput_bulk_list, axis=1)
        logger.info("Metagenomics Flowcell id: {} - Bulk creating CentOutput objects".format(self.flowcell_id))
        CentOutput.objects.bulk_create(list(centoutput_insert_list.values))

        df.set_index("tax_id", inplace=True)

        df = df.append(barcode_df)
        df.reset_index(inplace=True)
        # append previous data frame to current
        bar_queryset = CentOutputBarcoded.objects.filter(output__flowcell__id=self.flowcell_id,
                                                         output__task__id=self.flowcell_job_id).values()

        if bar_queryset:
            logger.info("Metagenomics Flowcell id: {} - Previous CentOutput data found".format(self.flowcell_id))
            to_create_bar_df = pd.DataFrame(list(bar_queryset))
            to_create_bar_df["temp"] = "Y"
            df["temp"] = "N"
            to_create_bar_df = df.append(to_create_bar_df)
            # Keep duplicated lines, these are lines with values that need updating in the database
            to_update_bar_df = to_create_bar_df[to_create_bar_df.duplicated(subset=["tax_id", "barcode"],
                                                                            keep=False)]
            to_update_bar_df.set_index(["barcode", "tax_id"], inplace=True)
            to_update_bar_df["updated_sum_unique"] = to_update_bar_df.groupby(["barcode",
                                                                               "tax_id"])["sum_unique"].sum()
            to_update_bar_df["updated_num_matches"] = to_update_bar_df.groupby(["barcode",
                                                                                "tax_id"])["num_matches"].sum()
            # to_update_bar_df = to_update_bar_df[~to_update_bar_df.index.duplicated(keep="first")]
            to_update_bar_df = to_update_bar_df[to_update_bar_df.index.duplicated(keep="first")]
            to_update_bar_df.reset_index(inplace=True)
            to_update_bar_df.apply(update_bar_values, args=(self.flowcell_job_id, self.flowcell_id),
                                   axis=1)

            # ###### Update existing barcode output ######
            to_create_bar_df.drop_duplicates(subset=["tax_id", "barcode"], inplace=True, keep=False)
            to_create_bar_df = to_create_bar_df[to_create_bar_df["temp"] == "N"]
            bulk_insert_list_bar = to_create_bar_df.apply(bulk_create_list, args=(job_master, self.flowcell_id),
                                                          axis=1)
            CentOutputBarcoded.objects.bulk_create(list(bulk_insert_list_bar.values))
            logger.info("Metagenomics Flowcell id: {} - Bulk creating CentOutputBarcoded objects".format(
                self.flowcell_id)
            )
            df.drop(columns=["temp"], inplace=True)

        else:
            bulk_insert_list_bar = df.apply(bulk_create_list, args=(job_master, self.flowcell_id),
                                            axis=1)
            CentOutputBarcoded.objects.bulk_create(list(bulk_insert_list_bar.values))
            logger.info("Metagenomics Flowcell id: {} - Bulk creating CentOutputBarcoded objects".format(
                self.flowcell_id)
            )

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
        # TODO rewrite into pandas mapping function, low priority, or magic datafame magic
        for key, value in lineages_taxidlist_dict.items():
            # get the ranks of a list of taxID. Returns a dict with keys of the taxID in the list of taxIDS,
            # value is the rank, i.e root, kingdom etc.
            # TODO vectorise this sucka
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

        # lin_df.fillna("NA")
        logger.info("Metagenomics Flowcell id: {} - Determining Subspecies".format(self.flowcell_id))
        # create new additional subspecies column, from strains column which has subspecies above it in taxa level
        lin_df["subStrainSpecies"] = lin_df["strain"].map(subspecies_determine)
        # merge new subspecies column with existing column
        lin_df["subspecies"].fillna(lin_df["subStrainSpecies"], inplace=True)
        # delete the new subspecies column
        self.delete_series(["subStrainSpecies"], lin_df)
        # ##############  Storing snakey links and nodes ################
        logger.info("Metagenomics Flowcell id: {} - Starting Sankey Calculations".format(self.flowcell_id))
        lin_df.head(n=20)
        # Order the DataFrame descending hierarchically, from kingdom to species
        merge_lin_df = pd.merge(lin_df, df, how="outer", left_index=True, right_on="tax_id")

        merge_lin_df.set_index("tax_id", inplace=True)

        columns = ["num_matches", "sum_unique", "barcode", "superkingdom", "phylum", "class", "order", "family",
                   "genus", "species"]
        sankey_lin_df = merge_lin_df[columns]
        # Baqckfill the dataframe nas with the lowest common entry, if family is missing, and genus
        #  isn't, family NaN becomes the genus entry
        sankey_lin_df = sankey_lin_df.fillna(axis=1, method="bfill")
        # Add the number of reads of the lowest common ancestor of a row to the dataframe
        # initialise the series that will become our sankey pandas DataFrame

        source = pd.Series()
        target = pd.Series()
        value = pd.Series()
        tax_id = pd.Series()
        barcode = pd.Series()
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
                barcode = barcode.append(sankey_lin_df["barcode"])
                sankey_lin_df["target_tax_level"] = target_clade
                target_tax_level = target_tax_level.append(sankey_lin_df["target_tax_level"])
        # Create a DataFrame of links
        source_target_df = pd.concat([source, target, value, tax_id, target_tax_level, barcode], axis=1)
        # Rename the columns
        source_target_df.columns = ["source", "target", "value", "tax_id", "target_tax_level", "barcode"]

        source_target_df["flowcell"] = flowcell
        source_target_df["job_master"] = job_master
        logger.info("Metagenomics Flowcell id: {} - Barcodes are:".format(self.flowcell_id))
        logger.info(source_target_df["barcode"].unique())
        # TODO DO sorting magic before the database deposition
        # Taxids from the current results that are in this df
        # current_links_taxids = set(source_target_df["tax_id"])
        # TAxids that have links from this analysis already
        prev_link_df = pd.DataFrame(
            list(SankeyLinks.objects.filter(flowcell__id=self.flowcell_id,
                                            task__id=self.flowcell_job_id).values()))
        if prev_link_df.empty:
            to_create_sank_df = source_target_df
        else:
            tax_id_list = list(prev_link_df["tax_id"])
            barcodes_list = list(prev_link_df["barcode"])
            tupley_tuple = list(zip(tax_id_list, barcodes_list))
            prev_links_mask = source_target_df[['tax_id', 'barcode']].agg(tuple, 1).isin(tupley_tuple)
            to_create_sank_df = source_target_df[~prev_links_mask]
        # Not yet linked is the taxIDs that don't have any result in the database
        # Create a subset dataframe that contains only links for new taxIDS
        # Create  a new series that contains the task record for this in each row
        # to_create_df["job_master"] = job_master

        # ## ### Create the SankeyLinks object #### # One for each link
        to_create_sank_df.reset_index(inplace=True)
        to_create_sank_df = to_create_sank_df[to_create_sank_df["value"] > 5]
        to_create_sank_df = to_create_sank_df.drop_duplicates(subset=["tax_id", "target_tax_level"])

        # Create a subset df containing links that already have reuslts
        # # add a new series, broadcast the flowcell object down it
        # to_create_df["flowcell"] = flowcell
        # to_update_df["flowcell"] = flowcell
        to_update_sank_df = pd.DataFrame()
        if not prev_link_df.empty:
            to_update_sank_df = pd.merge(source_target_df, prev_link_df, how="inner", on=["tax_id", "barcode",
                                                                                          "target_tax_level"])

        if not to_create_sank_df.empty:
            # Perform the apply

            sankey_link_insert_list = to_create_sank_df.apply(sankey_bulk_insert_list, axis=1)

            # Bulk create the objects
            SankeyLinks.objects.bulk_create(list(sankey_link_insert_list.values))

        if not to_update_sank_df.empty:
            to_update_sank_df["updated_value"] = to_update_sank_df["value_x"] + to_update_sank_df["value_y"]
            # TODO slowest point in code
            logger.info("Metagenomics Flowcell id: {} - Updating Sankey Links - This will take a minute...".format(
                self.flowcell_id)
            )
            to_update_sank_df.apply(sankey_bulk_bar_insert, args=(self.flowcell_id,
                                                                  self.flowcell_job_id), axis=1)
            logger.info("Metagenomics Flowcell id: {} - Finished updating sankey links".format(self.flowcell_id))

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
            lineages_insert_list = new_lineages_for_saving_df.apply(lineages_bulk_insert_list,
                                                                    axis=1)
            # Bulk create
            LineageValues.objects.bulk_create(list(lineages_insert_list.values))

        # Update the jobmaster object fields that are relevant
        # job_master.last_read = last_read_id
        # job_master.read_count = doc_no
        # job_master.save()
        metadata.reads_classified = 8008135
        metadata.save()
        # if iteration_count < 3:
        #     logger.info(iteration_count)
        #     timmy.sleep(60)
        # job_master.complete = True
        # job_master.save()
        # set running to false. This means client stops querying
        # MetaGenomicsMeta.objects.filter(meta_id=self.foreign_key).update(running=False)
        logger.info("Metagenomics Flowcell id: {} - Finished all reads in database".format(self.flowcell_id))

        start_time = metadata.timestamp.replace(tzinfo=None)
        end_time = timezone.now().replace(tzinfo=None)
        metadata.finish_time = str(end_time - start_time)
        metadata.save()
        # Update the jobmaster object fields that are relevant
        job_master.running = False
        job_master.last_read = self.last_read + self.chunk
        job_master.read_count = doc_no
        job_master.save()
        logger.info("Metagenomics Flowcell id: {} - Finished!".format(self.flowcell_id))
        logger.info("Metagenomics Flowcell id: {} - New last_read_id is - {}".format(self.flowcell_id, job_master.last_read))
        logger.info("Metagenomics Flowcell id: {} - Total CentOut lines - {}".format(self.flowcell_id, total_centout))
