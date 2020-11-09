"""Centrifuge.py
"""

import os
import subprocess
import sys
from collections import defaultdict
from io import StringIO

import numpy as np
import pandas as pd
from celery.utils.log import get_task_logger
from django.db.models import Sum
from django.utils import timezone
from ete3 import NCBITaxa

from metagenomics.centrifuge_validation import (
    separate_target_cent_output,
    create_mapping_result_objects,
    fetch_concat_reference_path,
    map_target_reads,
    save_mapping_results,
)
from metagenomics.models import (
    CentrifugeOutput,
    Metadata,
    DonutData,
)
from metagenomics.utils import (
    convert_species_to_subspecies,
    divd,
    create_centrifuge_models,
    delete_series,
    create_donut_data_models,
)
from minotourapp.celery import app
from minotourapp.utils import get_env_variable
from reads.models import JobMaster
from reads.utils import get_fastq_df

pd.options.mode.chained_assignment = None
logger = get_task_logger(__name__)


def run_centrifuge(flowcell_job_id, streamed_reads=None):
    """
    Run the metagenomics subprocess command, returning the data from it as a DataFrame.
    Parameters
    ----------
    flowcell_job_id: int
        The primary key of the flowcell ID
    streamed_reads: list of dict
        A list of dictionaries containing read information

    Returns
    -------
    pd.core.frame.DataFrame, int, int, int, pandas.core.frame.DataFrame, int, int
        Dataframe of metagenomics results, total output lines from metagenomics, last read primary key,
        total count of reads analysed, dataframe of any reads that identified as targets,
         number of reads with classifications, number of reads without classifications

    """
    # The JobMaster object
    task = JobMaster.objects.get(pk=flowcell_job_id)
    # The flowcell the reads are from
    flowcell = task.flowcell
    avg_read_length = int(flowcell.average_read_length)
    if avg_read_length == 0:
        logger.error(
            f"Average read length is zero Defaulting to 450, but this is an error."
        )
        avg_read_length = 1000
    if not streamed_reads and not isinstance(streamed_reads, list):
        read_count, last_read, fasta_df_barcode = get_fastq_df(
            flowcell_pk=int(flowcell.id),
            desired_yield=50,
            avg_read_len=avg_read_length,
            task=task,
        )
    else:
        last_read = task.last_read
        fasta_df_barcode = pd.DataFrame(streamed_reads)
        if not fasta_df_barcode.empty:
            fasta_df_barcode = fasta_df_barcode.rename(
                columns={"type": "read_type_id", "barcode": "barcode_id"}
            )
            fasta_df_barcode["type__name"] = fasta_df_barcode["read_type_id"]
        read_count = fasta_df_barcode.shape[0]
    print(fasta_df_barcode)
    if fasta_df_barcode.empty:
        return pd.DataFrame(), None, None, None, None, 0, 0
    print("Flowcell id: {} - number of reads found {}".format(flowcell.id, read_count))
    # Create a fastq string to pass to Centrifuge
    fasta_df_barcode["fasta"] = (
        ">read_id="
        + fasta_df_barcode["read_id"]
        + ",barcode="
        + fasta_df_barcode["barcode_name"]
        + "\n"
        + fasta_df_barcode["sequence"]
    )
    fastqs_data = "\n".join(list(fasta_df_barcode["fasta"]))
    print("Flowcell id: {} - Loading index and Centrifuging".format(flowcell.id))
    # Write the generated fastq file to stdin, passing it to the command
    # Use Popen to run the metagenomics command
    # The path to the metagenomics executable
    centrifuge_path = get_env_variable("MT_CENTRIFUGE")
    # The path to the Centrifuge Index
    index_path = get_env_variable("MT_CENTRIFUGE_INDEX")
    # The command to run metagenomics
    cmd = "perl " + centrifuge_path + " -f --mm -k 3 -p 3 -x " + index_path + " -"
    try:
        out, err = subprocess.Popen(
            cmd.split(),
            preexec_fn=lambda: os.nice(-10),
            stdout=subprocess.PIPE,
            stdin=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).communicate(input=str.encode(fastqs_data))
    except subprocess.SubprocessError as e:
        print(f"{e}, running with standard niceness index.")
        out, err = subprocess.Popen(
            cmd.split(),
            stdout=subprocess.PIPE,
            stdin=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).communicate(input=str.encode(fastqs_data))
    # The standard error
    if err:
        print("Flowcell id: {} - Centrifuge error!! {}".format(flowcell.id, err))
    # out is a bytestring so it needs decoding
    if not out:
        print(
            "Flowcell id: {} - No reads found or no metagenomics output."
            " Check above for error".format(flowcell.id)
        )
        task.running = False
        task.save()
        sys.exit(0)
    centrifuge_output = out.decode()
    # total number of lines of metagenomics output dealt with
    total_centrifuge_output = centrifuge_output.count("\n") - 1
    print(
        "Flowcell id: {} - number of metagenomics output lines is {}".format(
            flowcell.id, total_centrifuge_output
        )
    )
    # output fields is the column headers for the pandas data frame
    output_fields = ["readID", "seqID", "taxID", "numMatches"]
    # create the DataFrame from the output
    df = pd.read_csv(StringIO(centrifuge_output), sep="\t", usecols=output_fields)
    # split out the barcode_name from the readID string
    df = split_read_id_and_barcodes(df)
    individual_reads_classified = np.unique(df["readID"].values).size
    targets_df = separate_target_cent_output(df, task, fasta_df_barcode)
    # The number of reads we have any form of classification for
    reads_classified = np.unique(df[
        df["tax_id"].ne(0)
        ]["read_id"].values).size
    # The number of reads we have completely failed to classify
    reads_unclassified = np.unique(df[
        df["tax_id"].eq(0)
        ]["read_id"].values).size    # save the values
    # Get the metadata object. Contains the start time, end time and runtime of the task
    metadata, created = Metadata.objects.get_or_create(task=task)
    return df, individual_reads_classified, read_count, last_read, targets_df, reads_classified, reads_unclassified


def split_read_id_and_barcodes(df):
    """
    Split the combined read_id and barcode field into separate read_id and barcode columns
    Parameters
    ----------
    df: pd.core.frame.DataFrame
        Dataframe of metagenomics output results

    Returns
    -------
    pd.core.frame.DataFrame
        Dataframe of the output results with a split read id column
    """
    barcode_name_series = df["readID"].str.split(",").str[1].str.split("=").str[1]
    barcode_name_series = barcode_name_series.rename("barcode_name")
    df = pd.concat([df, barcode_name_series], axis=1)
    df = pd.concat(
        [df, df["readID"].str.split(",").str[0].str.split("=").str[1]], axis=1
    )
    # Add the column names
    df.columns = [
        "readID",
        "seqID",
        "tax_id",
        "numMatches",
        "barcode_name",
        "read_id",
    ]
    return df


def process_centrifuge_output(df, task):
    """
    Process the output of the metagenomics DataFrame, add in lineages of the classified taxa
    Parameters
    ----------
    df: pandas.core.frame.DataFrame
        Dataframe of the metagenomics output
    Returns
    -------
    pd.core.frame.DataFrame
        The results

    """
    barcodes = df["barcode_name"].unique().tolist()
    # Add an all reads string to the list
    barcodes.append("All reads")
    df["unique"] = np.where(df["numMatches"] == 1, 1, 0)  # condition  # True  # False
    # get np array of the taxids
    taxid_list = np.unique(df["tax_id"].values)
    # use the ncbi thing to get the species names for the tax ids in the tax ids list
    ncbi = NCBITaxa()
    taxid2name_df = pd.DataFrame.from_dict(
        ncbi.get_taxid_translator(taxid_list), orient="index", columns=["name"]
    )
    # insert the taxid2name dict items into the dataframe name columns
    df = pd.merge(df, taxid2name_df, how="outer", left_on="tax_id", right_index=True)
    # any NaNs replace with Unclassified
    df["name"].fillna("Unclassified", inplace=True)
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
    # Create a list of all the barcodes, one element for this list
    return df, task, barcodes


def calculate_barcoded_values(barcode_group_df, barcode_df, classified_per_barcode):
    """
    Return the main cent output data frame with the barcode only result concatenated on
    :param barcode_group_df: The dataframe of one barcode after a group by on all the results
    :param barcode_df: the main results data frame that will store the output of this function
    :param classified_per_barcode: The number of classified reads in a barcode
    :return: A new df with the result for barcode appended
    """
    # Get the barcodes for this run
    barcode = barcode_group_df["barcode_name"].unique()[0]
    # Get The total number of reads classified in this barcode from a dict containing all the values, keyed by barcode
    reads_classed_in_barcode = classified_per_barcode[barcode]
    # Group the dataframe by the tax id
    barcode_group_df_gb = barcode_group_df.groupby("tax_id")
    # calculate total unique matches for this species in this barcode as a series
    sum_unique = barcode_group_df_gb["unique"].sum()
    # calculate the total matches fot this species in this barcode as a series
    num_matches = barcode_group_df_gb.size()
    # Drop columns we don't use
    barcode_group_df.drop(
        columns=["readID", "seqID", "numMatches", "unique"], inplace=True
    )
    # Drop duplicates rows of species, as we have a summed number of matches form above
    barcode_group_df.drop_duplicates(
        subset=["barcode_name", "name"], keep="first", inplace=True
    )
    # Calculate the proportion
    barcode_group_df["proportion_of_classified"] = (
        barcode_group_df["num_matches"]
        .div(reads_classed_in_barcode)
        .mul(100)
        .round(decimals=3)
    )
    # Create a dataframe from the series of useful information we just calculated
    values_df = pd.concat(
        [
            sum_unique,
            num_matches,
            barcode_group_df["name"],
            barcode_group_df["proportion_of_classified"],
        ],
        axis=1,
    )
    barcode_df = barcode_df.append(values_df, sort=True)
    return barcode_df


def barcode_output_calculations(df, task):
    """
    Calculate the values (num_matches, sum_unique) grouped by barcode, and create one dataframe with these values
    Parameters
    ----------
    df: pd.core.frame.DataFrame
        The dataframe of results that we have
    task: reads.models.JobMaster
        The django ORM for this task
    Returns
    -------
    (pandas.core.frame.DataFrame, reads.models.JobMaster)
        Dataframe of all results including barcoded, The django ORM of the task

    """

    barcodes = np.unique(df["barcode_name"].values).tolist()
    # Add an all reads string to the list
    barcodes.append("All reads")
    classified_per_barcode = {
        barcode: CentrifugeOutput.objects.filter(
            task=task, barcode_name=barcode
        ).aggregate(Sum("num_matches"))["num_matches__sum"]
        for barcode in barcodes
    }
    if classified_per_barcode["All reads"] == None:
        # Set it to 0
        classified_per_barcode["All reads"] = 0
    gb_bc = df.groupby(["barcode_name"])
    # for the name and group_df in the barcode group by object
    # todo can be improved using counter
    for name, group in gb_bc:
        # If there is no value for that barcode
        if classified_per_barcode[name] is None:
            classified_per_barcode[name] = 0
        # set the value
        classified_per_barcode[name] += group.shape[0]
    classified_per_barcode["All reads"] += df.loc[~df.index.duplicated()].shape[0]
    barcode_df = pd.DataFrame()
    barcode_df = gb_bc.apply(
        calculate_barcoded_values, barcode_df, classified_per_barcode
    )
    # Reset the index
    barcode_df.reset_index(inplace=True)
    # rename the columns to useful names
    barcode_df.rename(columns={"unique": "sum_unique", 0: "num_matches"}, inplace=True)
    # Set the tax_id as the index
    barcode_df.set_index("tax_id", inplace=True)
    return barcode_df, task


def calculate_lineages_df(ncbi, df, tax_rank_filter, flowcell):
    """
    Insert any new lineages of species that are in our results into the database, and fetch all that we need for these
    results
    :param ncbi: The instance of the NCBITaxa class
    :param df: The data frame of our results
    :param tax_rank_filter: The order of our taxonomic ranks
    :param flowcell: The flowcell object
    :return lineages_Df: a DataFrame of all the lineages
    """
    # create new defaultDict which creates a default dict for missing values
    lineage_dict = defaultdict(lambda: defaultdict())
    # Get the current taxIDs in the Dataframe
    taxid_list = df["tax_id"].values
    # Check the taxIDs we have in the dataframe, to see if they have matching lineages already.
    # If not the set subtraction means that unique, lineage-less taxIDs end up in the not_prev_lineaged dict
    # Get the lineages for these taxIDs in the new metagenomics output,
    # returns a dict keyed by taxID,
    # value is listof taxIDS in lineage
    lineages_taxidlist_dict = ncbi.get_lineage_translator(taxid_list)
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
        lineage_dict[key] = {
            rank: taxid_species_lookup_dict[tax_id]
            for tax_id, rank in lineage_ranked.items()
            if rank in tax_rank_filter or key == tax_id
        }
    # create dataframe, the .T means that the dataframe is transposed, so dict key is the series title
    lineages_df = pd.DataFrame.from_records(lineage_dict).T
    # rename the no rank column to strain
    lineages_df = lineages_df.rename(columns={"no rank": "strain"})
    # if these series don't exist add them and populate with numpy NaN so there's no key errors
    for series in ["subspecies", "strain", "subStrainSpecies"]:
        if series not in lineages_df.keys():
            lineages_df[series] = np.NaN
    print("Flowcell id: {} - Determining Subspecies".format(flowcell.id))
    # create new additional subspecies column, from strains column which has subspecies above it in taxa level
    lineages_df["subStrainSpecies"] = lineages_df["strain"].map(
        convert_species_to_subspecies
    )
    # merge new subspecies column with existing column
    lineages_df["subspecies"].fillna(lineages_df["subStrainSpecies"], inplace=True)
    unclassified_row = {key: np.NaN for key in lineages_df.keys()}
    row = pd.DataFrame(unclassified_row, index=[0])
    lineages_df = lineages_df.append(row, sort=True)
    # delete the new subspecies column
    delete_series(["subStrainSpecies"], lineages_df)
    return lineages_df


def calculate_donut_data(df, flowcell, tax_rank_filter):
    """
    Calculate the donut data for the donut chart and table

    Parameters
    ----------
    df: pandas.core.frame.DataFrame
        The dataframe of the metagenomics output, aggregated
    flowcell: reads.models.Flowcell
        The flowcell ORM object
    tax_rank_filter: list
        List of taxonomic ranks

    Returns
    -------
    pd.core.frame.DataFrame
        data for the donut chart to save

    """
    # Merge the lineages dataframe onto the metagenomics output dataframe so each species has it's lineage
    print(df.keys())
    if "classy" in df.keys():
        df.rename(columns={"classy": "class"}, inplace=True)
    print(df.keys())
    print(df.head())
    data_df = df
    data_df.fillna("Unclassified", inplace=True)
    # Set the index to a multi index of all 7 taxonomic ranks
    data_df.set_index(tax_rank_filter, inplace=True)
    # group the dataframe by barcode
    gb_bc = data_df.groupby("barcode_name")
    # Initialise a dataframe to add onto
    donut_df = pd.DataFrame()
    # for the barcode_name and values (in dataframe form) for that barcode
    for name, group in gb_bc:
        # each taxonomic rank Example- "Kingdom", "Phylum" etc...
        for tax_rank in tax_rank_filter:
            # group the datframe by that level
            gb = group.groupby(level=tax_rank)
            # Get the number of matches for each member of that level
            series = gb["num_matches"].sum()
            # Create a dataframe of the matches for this tax rank
            temp_df = pd.DataFrame(series)
            # add a sum_unique column to the dataframe
            temp_df["sum_unique"] = gb["sum_unique"].sum()
            # set the tax rank in a column on that dataframe
            temp_df["tax_rank"] = tax_rank
            # Set the barcode for this set of data
            temp_df["barcode_name"] = name
            # Append the reuslts to the donut_df, so we have 7 sets (the tax rank) for each barcode
            donut_df = donut_df.append(temp_df)
    # As it only shows top ten, trim out any with less than 2 matches, as they definitely won't be making it
    donut_df = donut_df[donut_df["num_matches"] > 2]
    # name the index as name
    donut_df.index.name = "name"
    # Reset the index
    donut_df.reset_index(inplace=True)

    logger.info(
        "Flowcell id: {} - Bulk inserting new species donut data".format(flowcell.id)
    )
    donut_to_create_df = donut_df
    # Call output parser to insert the results dataframe
    return donut_to_create_df


def process_centrifuge_barcode_data(df, barcode_df, task, tax_rank_filter):
    """
    Process the metagenomics per barcode output dataframe append lineages for each species
    Parameters
    ----------
    df: pandas.core.frame.DataFrame
        The dataframe of All reads results
    barcode_df: pandas.core.frame.DataFrame
        Dataframe containing per barcode results
    task: reads.models.JobMaster
        The task ORM
    tax_rank_filter: list of str
        List of taxonomic ranks in hierarchical order

    Returns
    -------
    pd.core.frame.DataFrame
        The dataframe that we will be creating our output results from

    """
    flowcell = task.flowcell
    df.drop(
        columns=["readID", "seqID", "numMatches", "unique", "barcode_name", "read_id"],
        inplace=True,
    )

    # Make this the results for all reads in the metagenomics output
    df["barcode_name"] = "All reads"
    # Add the barcode dataframe onto the dataframe, so now we have all reads and barcodes
    df = df.append(barcode_df, sort=True)
    # drop all duplicated lines to keep only one for each entry, creating atomicity
    df.reset_index(inplace=True)
    df = df.set_index(["tax_id", "barcode_name"])
    df = df.loc[~df.index.duplicated()]
    # Set the flowcell object into a series, broadcast to each row
    df["flowcell"] = flowcell

    # same with task from Jobmaster
    # Get the task again so we don't over write the paused status
    task = JobMaster.objects.get(pk=task.id)
    df["task"] = task

    df.reset_index(inplace=True)
    # Create the lineages df, a df with the lineage of each species in the full results dataframe
    ncbi = NCBITaxa()
    lineages_df = calculate_lineages_df(ncbi, df, tax_rank_filter, flowcell)
    cent_to_create_df = pd.merge(
        df, lineages_df, how="inner", left_on="tax_id", right_index=True,
    )
    print("Merged in lineages")
    cent_to_create_df["task"] = task
    cent_to_create_df = cent_to_create_df[cent_to_create_df["barcode_name"] != "No"]
    return cent_to_create_df


def output_aggregator(task, new_data_df, is_donut):
    """
    Parse the metagenomics output or Donut data, update analyses values and return

    Parameters
    ----------
    task: reads.models.JobMaster
        The django ORM task
    new_data_df: pandas.core.frame.DataFrame
        The dataframe containing the parsed and aggregated metagenomics output results of this iteration
    is_donut: bool
        If analyses is for donut chart data

    Returns
    -------
    (pandas.core.frame.DataFrame, dict, int)
        Updated dataframe, number of total matches per barcode, the latest iteration count

    """

    metagenomics_task = task
    parsed_data = (
        CentrifugeOutput.objects.filter(
            task=metagenomics_task, latest=task.iteration_count
        )
        if not is_donut
        else DonutData.objects.filter(
            task=metagenomics_task, latest=task.iteration_count
        )
    )
    # Create a dataframe of the previous values
    parsed_data_df = pd.DataFrame.from_records(parsed_data.values())
    # class throws a key error cause of poor choices I made, so rename it classy
    if "class" in new_data_df.keys():
        new_data_df.rename(columns={"class": "classy"}, inplace=True)
    # If we have previous data
    if not parsed_data_df.empty:
        # Add one to the iteration count
        iteration_count = task.iteration_count + 1
        # If this is for CentrifugeOutput
        if not is_donut:
            # Merge the lsat Iterations values with this Iterations values
            merged_data_df = pd.merge(
                parsed_data_df, new_data_df, on=["barcode_name", "tax_id"], how="outer",
            )
            df_values = [
                ("superkingdom_x", "superkingdom_y", "superkingdom"),
                ("phylum_x", "phylum_y", "phylum"),
                ("classy_x", "classy_y", "classy"),
                ("order_x", "order_y", "order"),
                ("family_x", "family_y", "family"),
                ("genus_x", "genus_y", "genus"),
                ("species_x", "species_y", "species"),
                ("name_x", "name_y", "name"),
            ]
            for x in df_values:
                merged_data_df[x[0]].fillna(merged_data_df[x[1]], inplace=True)
            rename_values = {x[0]: x[2] for x in df_values}
            merged_data_df.rename(columns=rename_values, inplace=True)
        # If this is for a donut
        else:
            # Merge the last Iterations values with this Iterations values
            merged_data_df = pd.merge(
                parsed_data_df,
                new_data_df,
                how="outer",
                on=["barcode_name", "name", "tax_rank"],
            )

        # Set the number of matches and sum_unique NaNs to 0, NaNs present where there were missing values
        values = {
            "num_matches_x": 0,
            "num_matches_y": 0,
            "sum_unique_x": 0,
            "sum_unique_y": 0,
        }
        # Fill those NaNs with 0 so we can combine the columns
        merged_data_df.fillna(value=values, inplace=True)
        # Set the task column to the task objects
        merged_data_df["task"] = metagenomics_task
        # Combine the columns to get the updated number of matches
        merged_data_df["num_matches"] = (
            merged_data_df["num_matches_x"] + merged_data_df["num_matches_y"]
        )
        # Combine the columns to get the updated sum_unique
        merged_data_df["sum_unique"] = (
            merged_data_df["sum_unique_x"] + merged_data_df["sum_unique_y"]
        )
        # Group by barcode
        gb_bc = merged_data_df.groupby("barcode_name")
        # The number of matches per barcode as a dictionary
        classed_per_barcode = gb_bc["num_matches"].sum().to_dict()
        to_save_df = merged_data_df
    # Else this is the first iteration as we have no previous data
    else:
        # Set iteration count to one
        iteration_count = 1
        # Group by barcode name
        gb_bc = new_data_df.groupby("barcode_name")
        # The number of matches per barcode as a dictionary
        classed_per_barcode = gb_bc["num_matches"].sum().to_dict()
        to_save_df = new_data_df

    to_save_df["proportion_of_classified"] = to_save_df.apply(
        divd, args=(classed_per_barcode,), axis=1
    )
    # Set latest to the updated iteration count
    to_save_df["latest"] = iteration_count
    to_save_df["task"] = metagenomics_task
    return to_save_df, classed_per_barcode, iteration_count


def save_analyses_output_in_db(
    to_save_df, donut, task, iteration_count, classed_per_barcode=0
):
    """
    Save the aggregated output of the analyses into the database
    Parameters
    ----------
    to_save_df: pandas.core.frame.DataFrame
        Dataframe containing aggregated results to save
    donut: bool
        Whether the data is for the donut charts or not
    task: reads.models.JobMaster
        The task django ORM object
    iteration_count: int
        The latest iteration number
    classed_per_barcode: dict
        Dictionary containing the number of total metagenomics matches in a barcode

    Returns
    -------

    """
    if not donut:
        to_save_df.rename(columns={"classy": "class"}, inplace=True)
        to_save_df["proportion_of_classified"].fillna("Unclassified", inplace=True)
        to_save_df.fillna("Unclassified", inplace=True)

        # Create a series of CentrifugeOutput objects, one for each row in the dataframe
        to_bulk_save_series = to_save_df.apply(
            create_centrifuge_models, args=(classed_per_barcode,), axis=1,
        )
        # Bulk create all those objects in the series
        CentrifugeOutput.objects.bulk_create(
            to_bulk_save_series.values.tolist(), batch_size=1000
        )
        # Delete the values from the old iteration for the Centrifuge Output data
        CentrifugeOutput.objects.filter(task=task, latest__lt=iteration_count).delete()
    # If this is for a Donut
    else:
        # Create a series of DonutData objects, one for each row in the dataframe
        to_bulk_save_series = to_save_df.apply(
            create_donut_data_models, args=(task,), axis=1
        )
        # Bulk create all those objects in the series
        DonutData.objects.bulk_create(
            to_bulk_save_series.values.tolist(), batch_size=1000
        )
        logger.info(" Created donut data ")
        # Set the iteration count on the job master
        task = JobMaster.objects.get(pk=task.id)
        task.iteration_count = iteration_count
        task.save()
        # Delete the values from the old iteration from the Donut Data
        DonutData.objects.filter(task=task, latest__lt=iteration_count).delete()


def update_metadata_and_task(
    last_read, read_count, flowcell_job_id, total_centrifuge_output, reads_classified, reads_unclassified
):
    """
    Update metadata and task after job has finished
    Parameters
    ----------
    metadata: metagenomics.models.Metadata
        Metadata about this run
    last_read: int
        The primary key of the last fastqread we processed
    read_count: int
        The read count that we have
    flowcell_job_id: int
        The primary key of the job_master record in the database
    total_centrifuge_output: int
        The total metagenomics output
    reads_classified: int
        Number of reads with classifications
    reads_unclassified: int
        Number of reads without classifications

    Returns
    -------

    """
    task = JobMaster.objects.get(pk=flowcell_job_id)
    metadata, created = Metadata.objects.get_or_create(task=task)
    metadata.classified += reads_classified
    metadata.unclassified += reads_unclassified
    start_time = metadata.start_time.replace(tzinfo=None)
    # Get the end time
    end_time = timezone.now().replace(tzinfo=None)
    # Calculate the finish time
    metadata.finish_time = str(end_time - start_time)
    # Save the values
    metadata.save()
    # Refresh the task
    flowcell = task.flowcell
    # Set running to false
    task.running = False
    # OK - this returns the last record of all fastqs, not the last record of the chunk.
    task.last_read = last_read
    # Set the read count
    task.read_count = task.read_count + read_count
    # Save the new task values
    task.save()
    print(
        "Flowcell id: {} - New last_read_id is - {}".format(flowcell.id, task.last_read)
    )
    print(
        "Flowcell id: {} - Total CentOut lines - {}".format(
            flowcell.id, total_centrifuge_output
        )
    )
    print("Flowcell id: {} - Finished!".format(flowcell.id))


@app.task
def run_centrifuge_pipeline(flowcell_job_id, streamed_reads=None):
    """
    Run the metagenomics pipeline. TOdo could be split into a chain if slow
    Parameters
    ----------
    flowcell_job_id: int
        The primary key of the JobMaster record in the database
    streamed_reads: list of dict
        List of dictionaries pertaining to reads
    Returns
    -------

    """
    tax_rank_filter = [
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
    task = JobMaster.objects.get(pk=flowcell_job_id)
    df, total_centrifuge_output, read_count, last_read, targets_df, reads_classified, reads_unclassified = run_centrifuge(
        flowcell_job_id, streamed_reads
    )
    if df.empty:
        logger.info("No Centrifuge output, skipping iteration...")
        return
    df, task, barcodes = process_centrifuge_output(df, task)
    barcode_df, task = barcode_output_calculations(df, task)
    df = process_centrifuge_barcode_data(df, barcode_df, task, tax_rank_filter)
    to_save_df, classified_per_barcode, iteration_count = output_aggregator(
        task, df, False
    )
    save_analyses_output_in_db(
        to_save_df, False, task, iteration_count, classified_per_barcode
    )
    # ## mapping section
    target_region_df = create_mapping_result_objects(barcodes, task)
    path_to_reference = fetch_concat_reference_path(task)
    if not targets_df.empty:
        map_output_df = map_target_reads(
            task, path_to_reference, targets_df, to_save_df, target_region_df
        )
        if not map_output_df.empty:
            map_output_df.apply(save_mapping_results, axis=1, args=(task,))

    # #####  donut chart calculations
    donut_df = calculate_donut_data(df, task.flowcell, tax_rank_filter)
    to_save_df, classified_per_barcode, iteration_count = output_aggregator(
        task, donut_df, True
    )
    save_analyses_output_in_db(to_save_df, True, task, iteration_count)
    update_metadata_and_task(
        last_read, read_count, flowcell_job_id, total_centrifuge_output, reads_classified, reads_unclassified
    )
