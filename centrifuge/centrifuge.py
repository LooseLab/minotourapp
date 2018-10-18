"""
    Centrifuge.py
"""
import subprocess
from collections import defaultdict
from datetime import datetime
from io import StringIO

import numpy as np
import pandas as pd
from celery.utils.log import get_task_logger
from django.db.models import ObjectDoesNotExist
from django.utils import timezone
from ete3 import NCBITaxa

from centrifuge.models import CentOutput, LineageValues, MetaGenomicsMeta, SankeyLinks, CentOutputBarcoded
from jobs.models import JobMaster
from minotourapp.utils import get_env_variable
from reads.models import FastqRead

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
    Create barcoded objects and return as an element in a new series
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
    Update existing CentOutputBarcoded objects in the database
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
    apply to the dataframe and create sankeylinks objects for each row, return into a series
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
    Update the barcode objects in the database
    :param row: The dataframe row
    :param flowcell_id: The primary key for the flowcell the task is being run on
    :param flowcell_job_id: The Pk of the task ID
    :return: Nothing
    """
    SankeyLinks.objects.filter(flowcell__id=flowcell_id, task__id=flowcell_job_id,
                               target_tax_level=row["target_tax_level"],
                               barcode=row["barcode"], tax_id=row["tax_id"]).update(value=row["updated_value"])


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


def insert_new_lineages(ncbi, df, tax_rank_filter, flowcell):
    """
    Insert any new lineages of species that are in our results into the database, and fetch all that we need for these
    results
    :param ncbi: The instance of the NCBITaxa class
    :param df: The dataframe of our results
    :param tax_rank_filter: The order of our taxonomic ranks
    :param flowcell: The flowcell object
    :return:
    """
    # create new defaultDict which creates a default dict for missing values
    new = defaultdict(lambda: defaultdict())
    # Get the current taxIDs in the Dataframe
    taxid_list2 = df["tax_id"].values
    # Get the taxIDs of lineages that are in the database, so we can use existing lineages or create new ones
    already_lineaged_tax_ids = LineageValues.objects.values_list("tax_id", flat=True)

    # Check the taxIDs we have in the dataframe, to see if they have matching lineages already.
    # If not the set subtraction means that unique, lineage-less taxIDs end up in the not_prev_lineaged dict
    not_prev_lineaged = list(set(taxid_list2) - set(already_lineaged_tax_ids))

    # Get the lineages for these taxIDs in the new centrifuge output,
    # returns a dict keyed by taxID,
    # value is listof taxIDS in lineage
    lineages_taxidlist_dict = ncbi.get_lineage_translator(taxid_list2)
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
    logger.info("Flowcell id: {} - Determining Subspecies".format(flowcell.id))
    # create new additional subspecies column, from strains column which has subspecies above it in taxa level
    lin_df["subStrainSpecies"] = lin_df["strain"].map(subspecies_determine)
    # merge new subspecies column with existing column
    lin_df["subspecies"].fillna(lin_df["subStrainSpecies"], inplace=True)

    # delete the new subspecies column
    delete_series(["subStrainSpecies"], lin_df)
    return lin_df


def calculate_insert_sankey(lin_df, df, flowcell, tax_rank_filter, job_master):
    """

    :param lin_df:
    :param df:
    :param flowcell:
    :param tax_rank_filter:
    :param job_master:
    :return:
    """
    # ##############  Storing snakey links and nodes ################

    logger.info("Flowcell id: {} - Starting Sankey Calculations".format(flowcell.id))
    # pd.merge outer join to combine the rows on tax_id, one row for each tax_id
    merge_lin_df = pd.merge(lin_df, df, how="outer", left_index=True, right_on="tax_id")
    # set index to tax_id
    merge_lin_df.set_index("tax_id", inplace=True)

    # The columns we actually need
    columns = ["num_matches", "sum_unique", "barcode", "superkingdom", "phylum", "class", "order", "family",
               "genus", "species"]

    # Subset the dataframe into those columns
    sankey_lin_df = merge_lin_df[columns]

    # Backfill the dataframe nas with the lowest common entry, if family is missing, and genus
    #  isn't, family NaN becomes the genus entry
    sankey_lin_df = sankey_lin_df.fillna(axis=1, method="bfill")

    # initialise the series that will become our sankey pandas DataFrame
    source = pd.Series()
    target = pd.Series()
    value = pd.Series()
    tax_id = pd.Series()
    barcode = pd.Series()
    target_tax_level = pd.Series()
    # Iterate across the dataframe taxonnomic rank columns, to create the source, target, num reads series
    for index, rank in enumerate(tax_rank_filter):
        # The taxonomic rank of the source of this link
        source_trank = rank
        # The next index up, so for kingdom it would be phylum
        k = index + 1
        if k < len(tax_rank_filter):
            # The target taxonomic rank
            target_trank = tax_rank_filter[k]
            # Create a source series, append it to the bottom of the existing source series
            source = source.append(sankey_lin_df[source_trank])
            # Create a target series, append it to the bottom of the existing target series
            target = target.append(sankey_lin_df[target_trank])
            # Create a values series, append it to the bottom of the existing value series,
            # contains the number of matches
            value = value.append(sankey_lin_df["num_matches"])
            # Create a tax_id series, append it to the bottom of the existing tax_id series
            tax_id = tax_id.append(pd.Series(sankey_lin_df.index.values, index=sankey_lin_df.index))
            # Create a barcode series, append it to the bottom of the existing barcode series
            barcode = barcode.append(sankey_lin_df["barcode"])
            # Create a target tax level series, by broadcasting the value down the dataframe,
            # add it to the bottom of the existing series
            # Contains the taxonomic rank of the links target, useful for ordering later
            sankey_lin_df["target_tax_level"] = target_trank
            target_tax_level = target_tax_level.append(sankey_lin_df["target_tax_level"])

    # Create a DataFrame of links
    source_target_df = pd.concat([source, target, value, tax_id, target_tax_level, barcode], axis=1)

    # Rename the columns
    source_target_df.columns = ["source", "target", "value", "tax_id", "target_tax_level", "barcode"]

    # Create a series of flowcell objects, through broadcasting
    source_target_df["flowcell"] = flowcell
    # Create a series of Job_master objects through broadcasting
    source_target_df["job_master"] = job_master

    logger.info("Flowcell id: {} - Barcodes are:".format(flowcell.id))
    logger.info(source_target_df["barcode"].unique())

    # TODO DO sorting magic before the database deposition
    # Create a dataframe of the previous results for this metagenomics task
    prev_link_df = pd.DataFrame(
        list(SankeyLinks.objects.filter(flowcell__id=flowcell.id,
                                        task__id=job_master.id).values()))

    # If there no results, i.e this is the first iteration
    if prev_link_df.empty:
        # The entire dataframe needs to be added into the database
        to_create_sank_df = source_target_df
    # if there are results these need updating where we have the same link
    else:
        tax_id_list = list(prev_link_df["tax_id"])

        barcodes_list = list(prev_link_df["barcode"])
        # list of tuples
        tupley_tuple = list(zip(tax_id_list, barcodes_list))
        # Mask whether there is a link for this barcode and tax_id
        prev_links_mask = source_target_df[['tax_id', 'barcode']].agg(tuple, 1).isin(tupley_tuple)
        # subset where there is so these need dootdating
        to_create_sank_df = source_target_df[~prev_links_mask]

    # ## ### Create the SankeyLinks object #### # One for each link
    to_create_sank_df.reset_index(inplace=True)

    # Remove any in this iteration where the value is 4 or less, as these will not be in the top 30
    to_create_sank_df = to_create_sank_df[to_create_sank_df["value"] > 5]
    # Drop duplicates as we only need to create one
    to_create_sank_df = to_create_sank_df.drop_duplicates(subset=["tax_id", "target_tax_level"])

    # Create a second dataframe for values that need updating
    to_update_sank_df = pd.DataFrame()

    # If there are values that need updating
    if not prev_link_df.empty:
        # Merge inner so only rows that are in both dataframe are kept, these are the rows that need updating
        # The target taxa level makes sure the rows intersect to the correct rank
        to_update_sank_df = pd.merge(source_target_df, prev_link_df, how="inner", on=["tax_id", "barcode",
                                                                                      "target_tax_level"])

    # If there are values that need creating
    if not to_create_sank_df.empty:
        # Perform the apply to create the series to create
        sankey_link_insert_list = to_create_sank_df.apply(sankey_bulk_insert_list, axis=1)

        # Bulk create the objects
        SankeyLinks.objects.bulk_create(list(sankey_link_insert_list.values))

    # If there are values that need updating
    if not to_update_sank_df.empty:
        # Sum the two values from the merged previous and current dataframes together
        to_update_sank_df["updated_value"] = to_update_sank_df["value_x"] + to_update_sank_df["value_y"]
        logger.info("Flowcell id: {} - Updating Sankey Links - This will take a minute...".format(
            flowcell.id)
        )
        # TODO slowest point in code
        # Apply the update function to the dataframe
        to_update_sank_df.apply(sankey_bulk_bar_insert, args=(flowcell.id,
                                                              job_master.id), axis=1)
        logger.info("Flowcell id: {} - Finished updating sankey links".format(flowcell.id))


def run_centrifuge(flowcell_job_id):
    """

    Returns nothing.
    -------
     Run the centrifuge command, and open the subprocess to run it. Keep the process open and pipe the reads into
     it Write the reads incrementally to the processes stdin, capture the stdout and parse it into a useful format.
     Once parsed it's deposited into the database, updating any existing records. The organising key is a uuid,
     stored in self.foreign key.
    """

    job_master = JobMaster.objects.get(pk=flowcell_job_id)
    job_master.running = True
    job_master.save()

    flowcell = job_master.flowcell

    logger.info('Flowcell id: {} - Running centrifuge on flowcell {}'.format(flowcell.id, flowcell.name))
    logger.info('Flowcell id: {} - job_master_id {}'.format(flowcell.id, job_master.id))
    logger.info('Flowcell id: {} - last_read {}'.format(flowcell.id, job_master.last_read))
    logger.info('Flowcell id: {} - read_count {}'.format(flowcell.id, job_master.read_count))

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

    fastqs = FastqRead.objects.filter(run__flowcell=flowcell, id__gt=int(job_master.last_read)
                                      ).order_by('id')[:500]

    logger.info('Flowcell id: {} - number of reads found {}'.format(flowcell.id, fastqs.count()))
    try:
        metadata = MetaGenomicsMeta.objects.get(flowcell=flowcell, task=job_master)
    except (MetaGenomicsMeta.DoesNotExist, MetaGenomicsMeta.MultipleObjectsReturned):

        metadata = MetaGenomicsMeta(
            run_time=time,
            flowcell=flowcell,
            running=True,
            number_of_reads=99999,  # to be deleted?
            reads_classified=0,
            task=job_master
        )

        metadata.save()

    fastqs_list = fastqs.values_list('read_id', 'barcode_name', 'fastqreadextra__sequence')

    fastqs_data = "".join([str(">read_id=" + c[0] + ",barcode=" + c[1] + "\n" + c[2] + "\n") for c in fastqs_list
                           if c[2] is not None])

    # Write the generated fastq file to stdin, passing it to the command
    logger.info("Flowcell id: {} - Loading index and Centrifuging".format(flowcell.id))

    # Use Popen to run the centrifuge command
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE) \
        .communicate(input=str.encode(fastqs_data))

    # the number of centrifuge lines
    cent_out = out.decode()

    # total number of lines of centrifuge output dealt with
    total_centout += cent_out.count("\n") - 1

    logger.info("Flowcell id: {} - number of centrifuge output lines is {}".format(flowcell.id, total_centout))

    # output fields is the column headers for the pandas data frame
    output_fields = ["readID", "seqID", "taxID", "numMatches"]

    # create the DataFrame from the output
    df = pd.read_csv(StringIO(cent_out), sep="\t", usecols=output_fields)

    if df.empty:
        logger.info("Flowcell id: {} - No reads found or no centrifuge output".format(flowcell.id))
        return

    # create DataFrame column unique, where if the number of matches is 1, unique is 1, if matches is more than
    # 1 unique is 0, using numpy condition
    df["unique"] = np.where(
        df["numMatches"] == 1,  # condition
        1,  # True
        0  # False
    )

    df = pd.concat([df, df["readID"].str.split(",").str[1].str.split("=").str[1]], axis=1)

    df.columns = ["readID", "seqID", "tax_id", "numMatches", "unique", "barcode"]

    df = pd.concat([df, df["readID"].str.split(",").str[0].str.split("=").str[1]], axis=1)

    df.columns = ["readID", "seqID", "tax_id", "numMatches", "unique", "barcode", "read_id"]

    # get np array of the taxids
    taxid_list = np.unique(df["tax_id"].values)

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

    # name_df = df[df["name"].isin(temp_targets)]
    #
    # if not name_df.empty:
    #     map the target reads
    #     m = Metamap(name_df, flowcell.id, flowcell_job_id)
    #     m.map_the_reads()

    barcode_df = pd.DataFrame()

    # TODO vectorise these bad boys
    logger.info("Flowcell id: {} - The dataframe shape is {}".format(flowcell.id, df.shape))

    unique_barcode = set([c[1] for c in fastqs_list])

    if len(unique_barcode) > 1:

        logger.info("Flowcell id: {} - These are the barcodes in this set {}".format(flowcell.id, unique_barcode))

        barcode_df = gb_bc.apply(barcode_calculations, barcode_df)

        barcode_df.reset_index(inplace=True)

        barcode_df.rename(columns={"unique": "sum_unique", 0: "num_matches"}, inplace=True)

        barcode_df.set_index("tax_id", inplace=True)

    # delete these columns, no longer needed
    df.drop(columns=["readID", "seqID", "numMatches", "unique", "barcode", "read_id"], inplace=True)

    # remove duplicate rows in the data frame,so we only have one entry for each species
    df.drop_duplicates(keep="first", inplace=True)

    # filter out no ranks, like root and cellular organism
    tax_rank_filter = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

    # drop all duplicated lines to keep only one for each entry, Atomicity
    df.reset_index(inplace=True)

    df = df.drop_duplicates(keep="first")
    # Make this the results for all reads in the centrifuge output
    df["barcode"] = "All reads"
    # query the index get all the objects in the database for this analysis at this point
    queryset = CentOutput.objects.filter(task=job_master)
    # Get the previous tax ids in this analysis
    prev_df_tax_ids = list(queryset.values_list("tax_id", flat=True))

    # ###### BULK CREATE CENTOUTPUT OBJECTS ########
    # Give each row a flowcell object cell, as part of the flowcell_id series
    df["flowcell"] = flowcell

    # same with task from Jobmaster
    df["task"] = job_master

    # Split out the rows that have a tax_id that doesn't already have results
    cent_to_create = df[~df["tax_id"].isin(prev_df_tax_ids)]

    # apply row wise to append model representation objects into a series
    centoutput_insert_list = cent_to_create.apply(centoutput_bulk_list, axis=1)

    logger.info("Flowcell id: {} - Bulk creating CentOutput objects".format(flowcell.id))
    # Bulk create the objects
    CentOutput.objects.bulk_create(list(centoutput_insert_list.values))

    df.set_index("tax_id", inplace=True)
    # Add the barcode dataframe onto the dataframe, so now we have all reads and barcodes
    df = df.append(barcode_df)

    df.reset_index(inplace=True)
    # Get all the results for the barcoded entries, so each species will have multiple entries under different barcodes
    # TODO Split back into one table in database, barcoded results
    bar_queryset = CentOutputBarcoded.objects.filter(output__flowcell__id=flowcell.id,
                                                     output__task__id=flowcell_job_id).values()
    # If there are results
    if bar_queryset:
        logger.info("Flowcell id: {} - Previous CentOutput data found".format(flowcell.id))
        # Create a pandas dataframe
        to_create_bar_df = pd.DataFrame(list(bar_queryset))
        # These results will be removed from the dataframe
        to_create_bar_df["temp"] = "Y"
        # The original dataframe results will be kept
        df["temp"] = "N"
        # Append the two dataframes together
        to_create_bar_df = df.append(to_create_bar_df)
        # separate out lines with duplicated tax_id and barcodes, these are lines with values that need updating in the
        #  database, keep both lines
        to_update_bar_df = to_create_bar_df[to_create_bar_df.duplicated(subset=["tax_id", "barcode"],
                                                                        keep=False)]
        # set multi index to sort by tax_ids into each barcode
        to_update_bar_df.set_index(["barcode", "tax_id"], inplace=True)
        # Sum the sum unique for each species in each barcode
        to_update_bar_df["updated_sum_unique"] = to_update_bar_df.groupby(["barcode",
                                                                           "tax_id"])["sum_unique"].sum()
        # Sum the number of matches as above. Same as fr all reads calculations
        to_update_bar_df["updated_num_matches"] = to_update_bar_df.groupby(["barcode",
                                                                            "tax_id"])["num_matches"].sum()
        # Remove duplicated lines
        to_update_bar_df = to_update_bar_df[to_update_bar_df.index.duplicated(keep="first")]

        to_update_bar_df.reset_index(inplace=True)
        # ###### Update existing barcode output ######
        # Apply to update the existing entries in the databases
        to_update_bar_df.apply(update_bar_values, args=(flowcell_job_id, flowcell.id),
                               axis=1)
        # Drop duplicates of barcodes under each tax_id so we only create one entry to update
        to_create_bar_df.drop_duplicates(subset=["tax_id", "barcode"], inplace=True, keep=False)
        # Remove original df entries, so just barcoded data, not "All reads"
        to_create_bar_df = to_create_bar_df[to_create_bar_df["temp"] == "N"]
        # Create the model objects, returns a series
        bulk_insert_list_bar = to_create_bar_df.apply(bulk_create_list, args=(job_master, flowcell.id),
                                                      axis=1)
        # Bulk create the value
        CentOutputBarcoded.objects.bulk_create(list(bulk_insert_list_bar.values))
        logger.info("Flowcell id: {} - Bulk creating CentOutputBarcoded objects".format(
            flowcell.id)
        )
        # remove unnecessary temp series
        df.drop(columns=["temp"], inplace=True)

    else:
        # apply to create the models into a series
        bulk_insert_list_bar = df.apply(bulk_create_list, args=(job_master, flowcell.id),
                                        axis=1)
        # Bulk create
        CentOutputBarcoded.objects.bulk_create(list(bulk_insert_list_bar.values))

        logger.info("Flowcell id: {} - Bulk creating CentOutputBarcoded objects".format(
            flowcell.id)
        )

    lin_df = insert_new_lineages(ncbi, df, tax_rank_filter, flowcell)
    calculate_insert_sankey(lin_df, df, flowcell, tax_rank_filter, job_master)

    # ##################### Back to the non sankey calculations ####################
    # Update the metadata for the headers
    metadata.reads_classified = 8008135
    metadata.save()
    logger.info("Flowcell id: {} - Finished all reads in database".format(flowcell.id))

    # Get the start time #TODO remove the unused fields
    start_time = metadata.timestamp.replace(tzinfo=None)
    end_time = timezone.now().replace(tzinfo=None)
    metadata.finish_time = str(end_time - start_time)
    metadata.save()

    # Update the jobmaster object fields that are relevant
    job_master = JobMaster.objects.get(pk=flowcell_job_id)
    job_master.running = False

    if fastqs.count() > 0:

        job_master.last_read = fastqs[499].id

    job_master.read_count = job_master.read_count + fastqs.count()
    job_master.save()

    logger.info("Flowcell id: {} - New last_read_id is - {}".format(flowcell.id, job_master.last_read))
    logger.info("Flowcell id: {} - Total CentOut lines - {}".format(flowcell.id, total_centout))
    logger.info("Flowcell id: {} - Finished!".format(flowcell.id))
