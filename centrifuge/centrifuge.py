"""Centrifuge.py
"""

import os
import subprocess
from collections import defaultdict
from io import StringIO

import numpy as np
import pandas as pd
from celery.utils.log import get_task_logger
from django.db.models import ObjectDoesNotExist, Q
from django.db.models import Sum
from django.utils import timezone
from ete3 import NCBITaxa

from alignment.tasks_alignment import align_reads
from centrifuge.models import (
    CentrifugeOutput,
    LineageValue,
    Metadata,
    MappingResult,
    MappingTarget,
    DonutData,
)
from minotourapp.utils import get_env_variable
from reads.models import FastqRead, JobType, JobMaster
from reference.models import ReferenceInfo

pd.options.mode.chained_assignment = None
logger = get_task_logger(__name__)


def call_fetch_reads_cent(runs, chunk_size, last_read):
    """
    Call the fetch reads function to create a fastq for the process
    :param runs: List of all runs on the flowcell
    :param chunk_size: The target number of reads we want
    :param last_read: The previous last read we took
    :return:
    """
    # a list of the fastqs object to pass into the mapping function
    fasta_objects = list()
    # Initialise the data frame
    fastq_df_barcode = pd.DataFrame()
    while True:
        # Call fetch_reads_cent to actually query the database
        reads, last_read, read_count, fastasmchunk = fetch_reads_cent(
            runs, chunk_size, last_read
        )
        # Add fasta_objects chunk to the list
        fasta_objects += fastasmchunk
        # Append the reads_df to the fastq_df
        fastq_df_barcode = fastq_df_barcode.append(reads)
        # If we have enough reads (more than chunk_size) or no reads
        if len(fastq_df_barcode) >= chunk_size or len(reads) == 0:
            break
    # Update the read count with the number of reads we just fetched
    read_count = len(fastq_df_barcode)
    return fastq_df_barcode, last_read, read_count, fasta_objects


def fetch_reads_cent(runs, chunk_size, last_read):
    """
    Query the database for runs from a flowcell
    :param runs: The list of runs objects for this flowcell
    :param chunk_size: The target number of reads to have pulled back
    :param last_read: The id of the last read
    :return:
    """
    # Store the first read id and which run it's from this dict
    countsdict = dict()
    # The list of fastqread objects
    fasta_objects = list()
    # loop through the runs
    for run in runs:
        # fastqs = FastqRead.objects.filter(run=run, id__gt=int(last_read)).first()
        # Query the database for the reads objects
        fastqs = (
            FastqRead.objects.values_list("id")
            .filter(run=run, id__gt=int(last_read))
            .first()
        )
        # If there are fastqs store the object
        if fastqs:
            # Store the first id you get and the urn it's from
            countsdict[fastqs[0]] = run
    # initialise the count
    count = 1
    # initialise the a dataframe to append reads to
    fastq_df_barcode = pd.DataFrame()
    # If you have reads from more than one run
    if len(countsdict) > 1:
        # Start with first read entry from this run
        for entry in sorted(countsdict):
            # If we haven't done more iterations than we have dict entries
            if count < len(countsdict):
                # Do the query
                fastqs = FastqRead.objects.filter(
                    run=countsdict[entry],
                    id__gt=int(last_read),
                    id__lt=int(sorted(countsdict)[count]),
                )[:chunk_size]
                # Append new results to barcode
                fastq_df_barcode = fastq_df_barcode.append(
                    pd.DataFrame.from_records(
                        fastqs.values(
                            "read_id", "barcode_name", "sequence", "id"
                        )
                    )
                )

                # add to list of the fastQRead objects new objects
                fasta_objects += list(fastqs)
                # Set the new last read id
                last_read = fastq_df_barcode.id.max()
                # If we have enough reads, stop
                if len(fastq_df_barcode) >= chunk_size:
                    break
            count += 1
    # If we only have one run
    elif len(countsdict) == 1:
        """TODO This is risky and it breaks the logic - we may end up skipping reads from other runs"""
        # The key is the first fastqreads object retrieved
        mykey = list(countsdict.keys())[0]
        # Get the fastq objects from the database
        fastqs = FastqRead.objects.filter(
            run=countsdict[mykey], id__gt=int(last_read)
        )[:chunk_size]
        # append the results to the fastq dataframe
        fastq_df_barcode = fastq_df_barcode.append(
            pd.DataFrame.from_records(
                fastqs.values("read_id", "barcode_name", "sequence", "id")
            )
        )
        # Create a list of the fastQRead objects we will use in this objects
        fasta_objects += list(fastqs)
        # get the new last read id
        last_read = fastq_df_barcode.id.max()
    # get the read count of the new reads pulled out to later update the job_master total
    read_count = len(fastq_df_barcode)
    # Return everything to be used above
    return fastq_df_barcode, last_read, read_count, fasta_objects


# TODO del unused df and things


def calculate_barcoded_values(
    barcode_group_df, barcode_df, classified_per_barcode
):
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


def convert_species_to_subspecies(name):
    """
    Centrifuge doesn't output subspecies as a rank. Checks if subsp. is in a species name, and creates the
    corresponding rank and subspecies name
    Parameters
    ----------
    name - (type string) The name of the strain

    Returns - (type String) The subspecies as a string
    -------

    """
    # if the name is a string, not a float
    if type(name) is str:
        # Fairly obvious
        if "subsp." in name:
            basestring = name.split("subsp.", 1)[0]
            string = name.split("subsp.", 1)[1]
            string = string.split(" ")[1]
            subspecies = basestring + " subsp. " + string
            return subspecies
    else:
        pass


def create_lineage_models(row):
    """
    Apply function to the dataframe to populate a list with a LineageValue for each row
    :param row: The row of the dataframe
    :return: A LineageValue that is appended to a series
    """
    return LineageValue(
        superkingdom=row["superkingdom"],
        phylum=row["phylum"],
        tax_id=row["tax_id"],
        classy=row["class"],
        order=row["order"],
        family=row["family"],
        genus=row["genus"],
        species=row["species"],
        subspecies=row["subspecies"],
        strain=row["strain"],
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
    :param df: The data frame of our results
    :param tax_rank_filter: The order of our taxonomic ranks
    :param flowcell: The flowcell object
    :return lineages_Df: a DataFrame of all the lineages
    """
    # create new defaultDict which creates a default dict for missing values
    lineage_dict = defaultdict(lambda: defaultdict())
    # Get the current taxIDs in the Dataframe
    taxid_list = df["tax_id"].values
    # Get the taxIDs of lineages that are in the database, so we can use existing lineages or create new ones
    already_lineaged_tax_ids = LineageValue.objects.values_list(
        "tax_id", flat=True
    )

    # Check the taxIDs we have in the dataframe, to see if they have matching lineages already.
    # If not the set subtraction means that unique, lineage-less taxIDs end up in the not_prev_lineaged dict
    not_prev_lineaged = list(set(taxid_list) - set(already_lineaged_tax_ids))

    # Get the lineages for these taxIDs in the new centrifuge output,
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

    # some strains are strains of a subspecies. But ncbi only returns the terminal (leaf) which is the
    # strain and not the subspecies above it so this function splits out the subspecies and adds it to a new
    # series
    new_lineages_for_saving_df = lineages_df[
        lineages_df.index.isin(not_prev_lineaged)
    ]
    # if the df isn't empty
    if not new_lineages_for_saving_df.empty:
        # iterate over lineage dataframe to save to database #TODO rewrite into using mysql
        # set the index name so upon reset the added series has a name
        new_lineages_for_saving_df.index.name = "tax_id"

        new_lineages_for_saving_df.reset_index(inplace=True)
        # list to contain the bulk insert object
        lineages_insert_list = new_lineages_for_saving_df.apply(
            create_lineage_models, axis=1
        )
        # Bulk create
        LineageValue.objects.bulk_create(list(lineages_insert_list.values))

    logger.info("Flowcell id: {} - Determining Subspecies".format(flowcell.id))
    # create new additional subspecies column, from strains column which has subspecies above it in taxa level
    lineages_df["subStrainSpecies"] = lineages_df["strain"].map(
        convert_species_to_subspecies
    )
    # merge new subspecies column with existing column
    lineages_df["subspecies"].fillna(
        lineages_df["subStrainSpecies"], inplace=True
    )

    unclassified_row = {key: np.NaN for key in lineages_df.keys()}
    row = pd.DataFrame(unclassified_row, index=[0])
    lineages_df = lineages_df.append(row, sort=True)

    # delete the new subspecies column
    delete_series(["subStrainSpecies"], lineages_df)

    return lineages_df


def update_mapped_non_dangerous(row, task):
    """
    Update number of mapped reads that don't map to danger but do map to reference
    :param row: DataFrame row
    :param task: The related task for this analysis
    :return:
    """
    # Get the MappingResult object for the target species
    mapped_result = MappingResult.objects.get(
        task=task, species=row["name"], barcode_name=row["barcode_name"]
    )
    # Add on the new num_mapped
    mapped_result.num_mapped += row["num_mapped"]
    # Get the number of red reads to update the red proportion of classified
    rr = mapped_result.red_reads
    # Get that new value
    nm = mapped_result.num_mapped
    # Calculate the number of mapped against the number of matches
    mapped_result.mapped_proportion_of_classified = round(
        (nm / mapped_result.num_matches) * 100, 2
    )
    # Update the proportion of mapped red reads
    mapped_result.red_reads_proportion_of_classified = round(
        (rr / nm) * 100, 2
    )
    # Save the new values
    mapped_result.save()


def plasmid_mapping(row, species, fastq_list, flowcell, read_ids, fasta_df):
    """Map the reads for groups that have danger regions on plasmids

    :param pd.Series row: The row of the target dataframe for species with plasmid danger regions
    :param str species: The species name
    :param list fastq_list: The reads sequences
    :param reads.models.Flowcell flowcell: The flowcell for logging by flowcell_id
    :param read_ids: reads ids of the reads being mapped
    :param fasta_df: dataframe containing the sequence data for the centrifuge classified reads that have filtered to
    target reads
    :return plasmid_map_df: A list of dicts containing the information about the plasmid mappings
    """
    read_ids = read_ids

    flowcell = flowcell

    species = species.replace(" ", "_")

    reference_name = species + "_" + row["name"]

    logger.info(reference_name)

    reference = ReferenceInfo.objects.filter(name=reference_name).first()
    logger.info(f"refernce not first is {ReferenceInfo.objects.filter(name=reference_name)}\n")

    logger.info(f"reference is {reference}\n")
    # Get the jobtype for this target species mapping, so we can find the mapping job
    map_job_type = JobType.objects.get(name="Other")
    # Create or get the mapping job for this target species mapping
    mapping_task, created = JobMaster.objects.get_or_create(
        job_type=map_job_type, flowcell=flowcell, reference=reference
    )
    # Call the align reads function from the web tasks-alignment file so we can map the reads against the reference
    align_reads(fastq_list, mapping_task.id, fasta_df)
    # Get the output of the align reads function
    map_output = PafStore.objects.filter(
        job_master=mapping_task, qsn__in=read_ids
    )

    logger.info(
        "Flowcell id: {} - Plasmid mapping map output is {}".format(
            flowcell.id, map_output
        )
    )
    # The column names for the dataframe created by the align reads output
    # If there is no mapping output from the align reads function
    if not map_output:
        logger.info(
            "Flowcell id: {} - No mappings for plasmid {} on species {}".format(
                flowcell.id, row["name"], species
            )
        )
        return np.NaN
    # If there is mapping output
    else:
        # Create a dataframe
        plasmid_map_df = pd.DataFrame(list(map_output.values()))
        # Set the column names
        plasmid_map_df.rename(
        columns={
            "qsn": "query_seq_name",
            "qsl": "query_seq_len",
            "qs": "query_start",
            "qe": "query_end",
            "rs": "rel_strand",
            "tsn_id": "target_seq_name",
            "tsl": "target_seq_length",
            "ts": "target_start",
            "te": "target_end",
            "nrm": "num_residue_matches",
            "abl": "alignment_block_length",
            "mq": "mapping_qual",
        },
        inplace=True,
    )

        # TODO NOTE that atm, we're just counting any map that reads to a plasmid as a red read

        logger.info(
            "Flowcell id: {} - This many reads mapped evilly on this reference {} for species {}".format(
                flowcell.id, plasmid_map_df.shape[0], species
            )
        )
        # return our new plasmid dataframe, annoyingly each one becomes an element in a series
        return plasmid_map_df


def update_mapped_red_values(row, task):
    """
    Update the values that we have in the database for each species after getting results for this iteration
    :param row: The results_df row
    :param task: The task model object
    :return:
    """
    # Get the mapping result objects for this species
    mapped = MappingResult.objects.get(
        species=row["species"], task=task, barcode_name=row["barcode_name"]
    )
    # Get the number of mapped to this target species
    nm = mapped.num_mapped
    # Add on the new red reads values
    mapped.red_reads += row["red_num_matches"]
    # Get the proportion of classified (red reads against number of mapped)
    mapped.red_reads_proportion_of_classified = round(
        (mapped.red_reads / nm) * 100, 2
    )
    # Save the updated values
    mapped.save()


def map_all_the_groups(
    target_species_group_df,
    group_name,
    flowcell,
    gff3_df,
    targets_results_df,
    task,
    fastas,
):
    """
    Map the reads from the target data frames, after they've been grouped by species, so this is run once fpor each species
    :param target_species_group_df: A data frame that contains reads from only one species
    :param group_name: The name of the species in this group
    :param flowcell: The flowcell model from where these reads came from
    :param gff3_df: The data frame containing the target regions
    :param targets_results_df: The data frame containing the new targets mapping against reference
    :param task: The task object
    :param fastas: A list of the fastqread objects that we used in this iteration
    :return red_df: Any reads that map to plasmids and their information as a dataframe
    """
    logger.info("Target species df is :")
    logger.info(target_species_group_df)
    logger.info(
        "Flowcell id: {} - species being mapped is {}".format(
            flowcell.id, group_name
        )
    )
    # Get the target region for this species
    species_regions_df = gff3_df[gff3_df["species"] == group_name]
    # Get the species that is in this group
    species = group_name.replace(" ", "_")

    user_id = task.flowcell.owner

    try:
        # Get the reference for this species
        references = ReferenceInfo.objects.filter(
            Q(owner_id=user_id) | Q(private=False)
        ).filter(name=species)[0]
        logger.info("References for this metganomics")
        logger.info(references)
        logger.info("\n")
    except ObjectDoesNotExist:

        logger.info(
            "Flowcell id: {} - No reference found for species {}!".format(
                flowcell.id, species
            )
        )

        return

    # Get the read_ids for the reads that centrifuge has classified as target species
    read_ids = target_species_group_df["read_id"].values.tolist()

    # Filter the list of FastqRead objects so we only have the objects that have the desired read_ids
    fastqs = list(
        filter(lambda fastq: fastq.read_id in read_ids, list(fastas))
    )

    # Get the metagenomics mapping job type
    map_job_type = JobType.objects.get(name="Other")
    # Get or create the task for this mapping
    mapping_task, created = JobMaster.objects.get_or_create(
        job_type=map_job_type, flowcell=flowcell, reference=references
    )
    if created:
        logger.info(f"Created mapping target object for reference {references}")

    b = [
        {"read_id": fastq.read_id, "sequence": fastq.sequence}
        for fastq in fastqs
    ]

    fasta_df = pd.DataFrame(b)

    # Call the align_reads function to perform the mapping, funciton found in web tasks_alignment
    align_reads(fastqs, mapping_task.id, fasta_df)

    # Get the output for this mapping (if any) from the database
    map_output = PafStore.objects.filter(
        job_master=mapping_task, qsn__in=read_ids
    )

    logger.info(
        "Flowcell id: {} - minimap output {} ".format(flowcell.id, map_output)
    )

    # Initialise a dataframe to store the results of mapping to the plasmid
    plasmid_red_df = pd.DataFrame()

    # Get the target regions for plasmids for this species
    plasmid_df = species_regions_df[
        species_regions_df["gff_line_type"] == "plasmid"
    ]

    # If there are target regions for plasmids
    if not plasmid_df.empty:

        logger.info(
            "Flowcell id: {} - Mapping reads to plasmids for species {} ".format(
                flowcell.id, species
            )
        )
        # Map the reads against the target regions
        logger.info("plasmid df is")
        logger.info(plasmid_df)

        plasmid_red_series = plasmid_df.apply(
            plasmid_mapping,
            args=(species, fastas, flowcell, read_ids, fasta_df),
            axis=1,
        )
        # If a plasmid has no mappings, it results in a nan being returned into the series, so drop it so we only have
        # mapping results. IF no plasmids have mappings, the series becomes empty
        plasmid_red_series.dropna(inplace=True)

        logger.info(
            "Flowcell id: {} - plasmid mapping output is {}".format(
                flowcell.id, plasmid_red_series.head()
            )
        )
        # If there are mappings to a plasmid
        if not plasmid_red_series.empty:
            try:
                # for the index and the series (containg plasmid mapping dataframes as results)
                # of plasmid mapping results
                for ind, df_plas in plasmid_red_series.iteritems():
                    # error trapping, if for some reason we have an ndarray instead of a df (shouldn't)
                    if type(df_plas) == np.ndarray:
                        # make the ndarray into a dataframe
                        df_append = pd.DataFrame(df_plas)
                        # set the column names
                        # df_append.columns = ["alignment_block_length", "flowcell_id", "id", "job_master_id",
                        #                      "mapping_qual", "num_residue_matches",
                        #                      "query_end", "query_start", "query_seq_len", "read_id", "read_pk",
                        #                      "reference_id", "rel_strand",
                        #                      "target_end", "target_start", "target_seq_len", "target_seq_name"]

                        df_append.rename(
                            columns={
                                "qsn": "query_seq_name",
                                "qsl": "query_seq_len",
                                "qs": "query_start",
                                "qe": "query_end",
                                "rs": "rel_strand",
                                "tsn_id": "target_seq_name",
                                "tsl": "target_seq_length",
                                "ts": "target_start",
                                "te": "target_end",
                                "nrm": "num_residue_matches",
                                "abl": "alignment_block_length",
                                "mq": "mapping_qual",
                            },
                            inplace=True,
                        )
                        # append the plasmid mapping results to a consolidated all plasmid mappings dataframe
                        plasmid_red_df = plasmid_red_df.append(
                            df_append, sort=True
                        )
                    # if it is a dataframe (it should be)
                    else:
                        # append the plasmid mapping results to a consolidated all plasmid mappings dataframe
                        plasmid_red_df = plasmid_red_df.append(
                            df_plas, sort=True
                        )
                logger.info(
                    "Flowcell id: {} - plasmid mapping output is {}".format(
                        flowcell.id, plasmid_red_df.head()
                    )
                )
            # Catch and attribute error
            except AttributeError as e:
                logger.info(e)
                return

    # if there is mapping output (non plasmid) for the reference
    if map_output:
        # create a dataframe of the results
        map_df = pd.DataFrame(list(map_output.values()))
    # Else stop and return. There is no output if there are no mappings
    else:
        logger.info(
            "Flowcell id: {} - No minimap output (No reads mapped) for species {}".format(
                flowcell.id, species
            )
        )
        return

    # The column names
    # columns = ["alignment_block_length", "flowcell_id", "id", "job_master_id", "mapping_qual", "num_residue_matches",
    #            "query_end", "query_start", "query_seq_len", "read_id", "read_pk", "reference_id", "rel_strand",
    #            "target_end", "target_start", "target_seq_len", "target_seq_name"]
    # # set the column name to above
    # map_df.columns = columns

    map_df.rename(
        columns={
            "qsn": "query_seq_name",
            "qsl": "query_seq_len",
            "qs": "query_start",
            "qe": "query_end",
            "rs": "rel_strand",
            "tsn_id": "target_seq_name",
            "tsl": "target_seq_length",
            "ts": "target_start",
            "te": "target_end",
            "nrm": "num_residue_matches",
            "abl": "alignment_block_length",
            "mq": "mapping_qual",
        },
        inplace=True,
    )

    map_df["num_residue_matches"] = map_df["num_residue_matches"].astype(
        np.int64
    )
    # Filter out low quality matches or the
    map_df = map_df.query("mapping_qual >= 40 & num_residue_matches >= 200")

    # See whether the start or end of a mapping falls into the region
    boolean_df = species_regions_df.apply(
        falls_in_region, args=(map_df,), axis=1
    )
    # Get whether it's inside the boundaries or not
    boolean_df = boolean_df.any()
    # slice the dataframe to get the Red reads, reads that fall within boundaries
    red_df = map_df[boolean_df]
    logger.info(
        "Flowcell id: {} - The target df before plasmid mapping results is {}".format(
            flowcell.id, red_df
        )
    )
    # If there is output from the plasmid mapping, put it into the consolidated target reads containing dataframe
    if not plasmid_red_df.empty and (
        type(plasmid_red_df) != pd.core.series.Series
    ):
        # Filter out low mapping qualities
        plasmid_red_df = plasmid_red_df.query(
            "mapping_qual >= 40 & num_residue_matches >= 200"
        )
        red_df = red_df.append(plasmid_red_df, sort=True)

    logger.info(
        "Flowcell id: {} - The target df after plasmid mapping results is {}".format(
            flowcell.id, red_df
        )
    )

    logger.info(
        "Flowcell id: {} - This many reads mapped evilly on this reference {} for species {}".format(
            flowcell.id, red_df.shape[0], species
        )
    )

    logger.info(
        "Flowcell id: {} - This many reads mapped elsewhere on this reference {}".format(
            flowcell.id, map_df.shape[0]
        )
    )
    # Inner join Merge the mapped reads into the original dataframe, so we have only
    # the information for reads that mapped
    update_num_mapped_df = pd.merge(
        targets_results_df,
        map_df,
        how="inner",
        left_on="read_id",
        right_on="query_seq_name",
    )
    # Deep copy the dataframe so we can edit it without the changes being reflected in the original, one
    # to split by barcode, one to be all reads
    update_num_mapped_df_non_barcode = update_num_mapped_df.copy(deep=True)
    # Set the barcode to all reads as this dataframe isn't split by barcode
    update_num_mapped_df_non_barcode["barcode_name"] = "All reads"
    # Update the number mapped by getting the number of dataframe rows
    update_num_mapped_df_non_barcode[
        "num_mapped"
    ] = update_num_mapped_df_non_barcode.shape[0]
    # Drop duplicates on species
    update_num_mapped_df_non_barcode.drop_duplicates(
        subset="name", inplace=True
    )
    # Apply to update the number mapped for each mapping result in the database
    update_num_mapped_df_non_barcode.apply(
        update_mapped_non_dangerous, args=(task,), axis=1
    )
    # Update the barcoded num mapped but none dangerous reads
    update_num_mapped_df_barcoded = update_num_mapped_df
    # Group by barcode and species
    gb_bc_nm = update_num_mapped_df_barcoded.groupby(["barcode_name", "name"])
    # set index on the original so we can set the values from the group_by onto the dataframe
    update_num_mapped_df_barcoded.set_index(
        ["barcode_name", "name"], inplace=True
    )
    # set the number of mapped to the size of the group by
    update_num_mapped_df_barcoded["num_mapped"] = gb_bc_nm.size()
    # Reset the index
    update_num_mapped_df_barcoded.reset_index(inplace=True)
    # Drop duplicates on name and barcode, so only species unique in barcode
    update_num_mapped_df_barcoded.drop_duplicates(
        subset=["name", "barcode_name"], inplace=True
    )
    # Update the num mapped for the target species for each barcode
    update_num_mapped_df_barcoded.apply(
        update_mapped_non_dangerous, args=(task,), axis=1
    )

    logger.info(
        "Flowcell id: {} - Finished updating mapped reads".format(flowcell.id)
    )

    return red_df


def falls_in_region(row, map_df):
    """
    Does this reads mapping fall within any of the regions we have defined?
    :param row: The gff3 dataframe row
    :param map_df: The results of the mapping step in a pandas dataframe
    :return:
    """
    # TODO which region?
    # If the start is higher than the region start
    start_low = map_df["target_start"] > row["start"]
    # If the start is lower than the region end
    start_high = map_df["target_start"] < row["end"]
    # If the end is higher than the region end
    end_low = map_df["target_end"] > row["start"]
    # If the end is lower than the region end
    end_high = map_df["target_end"] < row["end"]
    # Concat the series together, create "truth" dataframe
    bool_df = pd.concat(
        [start_low, start_high, end_low, end_high], axis=1, ignore_index=True
    )
    # Keep where the start is in the region,
    bool_df["keep"] = np.where(
        (bool_df[0] & bool_df[1]) | (bool_df[2] & bool_df[3]), True, False
    )
    return bool_df["keep"]


def update_targets_no_mapping(row, task):
    """
    Update the targets to reflect the number of matches, if we've classified reads that don't map

    :param row: The target dataframe row
    :param task: The task database model object
    :return: None
    """
    # Get the Mapping result object for each target reads
    obj = MappingResult.objects.get(
        tax_id=row["tax_id"], task=task, barcode_name=row["barcode_name"]
    )
    # Set the number of matches
    obj.num_matches = row["num_matches"]
    # Set the number of unique matches
    obj.sum_unique = row["sum_unique"]
    # Get the number of mapped from teh database
    nm = obj.num_mapped
    # Get the proportion of classified
    obj.proportion_of_classified = row["prop_classed"]
    # Calculate the new proportion of mapped against total reads classified
    obj.mapped_proportion_of_classified = round(
        nm / row["num_matches"] * 100, 3
    )
    # Save the new values
    obj.save()


def calculate_num_matches_update(
    target_df, task, num_matches_targets_barcoded_df, classed_per_bar
):
    """
    Update the number of matches and sum_unique for target reads
    :param target_df: The dataframe containing the target reads
    :param task: The task model object
    :param num_matches_targets_barcoded_df: The number of matches per barcode for each target species in a dataframe
    :param classed_per_bar: The number of reads classified per barcode as a dict
    :return:
    """

    # create a dataframe for barcoded calculations
    barcode_df = target_df.reset_index()
    # create a dataframe for non barcoded calculations
    all_reads_df = target_df.copy(deep=True)

    all_reads_df.reset_index(inplace=True)
    # Set the barcode as all reads, as it is all reads irrespective of barcoding
    all_reads_df["barcode_name"] = "All reads"
    # Drop duplicates of species in the dataframe
    all_reads_df.drop_duplicates(subset="name", inplace=True)
    # If this dataframe is not empty, it means it isn't the first iteration
    if not num_matches_targets_barcoded_df.empty:
        # Merge the two dataframes on barcode and target idea, with a left join,
        # so all values are kept in left (target species rows) dataframe
        all_reads_df = pd.merge(
            all_reads_df,
            num_matches_targets_barcoded_df,
            how="left",
            on=["tax_id", "barcode_name"],
        )
        # set the num_matches_y to 0 for any columns that haven't got a match
        all_reads_df["num_matches_y"].fillna(0, inplace=True)
        # set the sum unqiue y to 0 for any columns that haven't got a match
        all_reads_df["sum_unique_y"].fillna(0, inplace=True)
        # Combine the two series to get a new number matches
        all_reads_df["num_matches"] = (
            all_reads_df["num_matches_x"] + all_reads_df["num_matches_y"]
        )
        # Combine the two series to get a new sum unique matches
        all_reads_df["sum_unique"] = (
            all_reads_df["sum_unique_x"] + all_reads_df["sum_unique_y"]
        )
    # Calculate the proportion of classifed for this species against all reads classified
    all_reads_df["prop_classed"] = (
        all_reads_df["num_matches"]
        .div(classed_per_bar["All reads"])
        .mul(100)
        .round(decimals=3)
    )

    # Apply the update_targets_no_mapping function to the dataframe to update the targets num_matches and sum_unique
    all_reads_df.apply(update_targets_no_mapping, args=(task,), axis=1)
    # Barcoded calculations now
    # Group the barcoded datframe by barcode_name and species name
    gb_bc = barcode_df.groupby(["barcode_name", "name"])
    # Set the index so we can map values back from the group_by back onto the dataframe by index
    barcode_df.set_index(["barcode_name", "name"], inplace=True)
    # Get the number of matches for each species in a barcode by the size of it's respective group
    barcode_df["num_matches"] = gb_bc.size()
    # Get the sum unique by the sum of unique series in each respective group
    barcode_df["sum_unique"] = gb_bc["unique"].sum()
    # Reset index
    barcode_df.reset_index(inplace=True)
    # If this isn't the first iteration, and we have previous data
    if not num_matches_targets_barcoded_df.empty:
        # merge the two dataframe together, keeping all the target rows
        barcode_df = pd.merge(
            barcode_df,
            num_matches_targets_barcoded_df,
            how="left",
            on=["tax_id", "barcode_name"],
        )
        # fill nans with 0
        barcode_df["num_matches_y"].fillna(0, inplace=True)
        # fill nans with 0
        barcode_df["sum_unique_y"].fillna(0, inplace=True)
        # Sum the two series to get the new number of matches value
        barcode_df["num_matches"] = (
            barcode_df["num_matches_x"] + barcode_df["num_matches_y"]
        )
        # Sum the two series to get the new sum unique matches value
        barcode_df["sum_unique"] = (
            barcode_df["sum_unique_x"] + barcode_df["sum_unique_y"]
        )
    # Set the total matches in barcode series, contains total matches in that barcode
    barcode_df["total_in_barcode"] = barcode_df["barcode_name"].map(
        classed_per_bar
    )
    # Calculate the proportion classed matches for this species against all matches in a barcode
    barcode_df["prop_classed"] = (
        barcode_df["num_matches"]
        .div(barcode_df["total_in_barcode"])
        .mul(100)
        .round(decimals=3)
    )
    # Drop duplicates of barcode_name and tax_id, so species are unique inside each barcode
    barcode_df.drop_duplicates(subset=["barcode_name", "tax_id"], inplace=True)
    # Update the target species objects num_matches and sum_unique fields
    barcode_df.apply(update_targets_no_mapping, args=(task,), axis=1)


def create_mapping_results_objects(task, flowcell, target_set):
    """
    Create Mapping result objects for each target species in the database, one set for each barcode
    :param task: The metagenomics task
    :param flowcell: The flowcell django queryset object, which the reads came from
    :param target_set: The target set of regions, used to find the Mapping target values from DB
    :return:
    """
    # Get the targets set in the gff files and uploaded
    target_regions_df = pd.DataFrame(
        list(MappingTarget.objects.filter(target_set=target_set).values())
    )
    # If there are no target regions found for the set name
    if target_regions_df.empty:
        logger.info(
            "Flowcell id: {} - No set of target regions found by that set name".format(
                flowcell.id
            )
        )
        return
    # Get np array of the target species
    # TODO split into function so code is DRYER
    target_species = target_regions_df["species"].unique()
    # Get np array of their tax_ids
    target_tax_id = target_regions_df["tax_id"].unique()
    # Combine into tuple, one for each
    target_tuple = list(zip(target_species, target_tax_id))
    # create one MappingResults entry for each target species
    for target in target_tuple:
        obj, created = MappingResult.objects.get_or_create(
            task=task,
            species=target[0],
            tax_id=target[1],
            barcode_name="All reads",
            defaults={"red_reads": 0, "num_mapped": 0},
        )
        if created:
            logger.info(
                "Flowcell id: {} - Result Mapping object created for {}".format(
                    flowcell.id, obj.species
                )
            )
    return target_regions_df


def map_the_reads(
    name_df,
    task,
    flowcell,
    num_matches_targets_barcoded_df,
    targets,
    class_per_bar,
    fastas,
):
    """
    Map the reads is called on the dataframe of targets that have been split out, splits them into a group by name, and
    applies the map_all_the_groups function, which returns any reads that map to target areas in a data frame.
    :param name_df: The targets dataframe
    :param task: The task model object
    :param flowcell: The flowcell model object
    :param num_matches_targets_barcoded_df: The number of matches per barcode for each target species in a data frame
    :param targets: The targets in this set of targets
    :param class_per_bar: The number of reads classified per barcode as a dict
    :param fastas: List of fastqread objects for this iteration
    :return:
    """
    logger.info(
        "Flowcell id: {} - Mapping the target reads".format(flowcell.id)
    )
    # get the targets dataframe
    targets_df = name_df
    # TODO currently hardcoded below
    target_set = task.target_set
    # if there are no targets identified by centrifuge in this iteration
    if targets_df.empty:
        logger.info(
            "Flowcell id: {} - No targets in this batch of reads".format(
                flowcell.id, targets_df.shape
            )
        )
        # Create the mapping targets
        create_mapping_results_objects(task, flowcell, target_set)

        return

    # Reset the index
    targets_df.reset_index(inplace=True)

    logger.info(
        "Flowcell id: {} - targets_df shape is {}".format(
            flowcell.id, targets_df.shape
        )
    )

    # Create the mapping targets, calls get or create so they are only created once
    target_regions_df = create_mapping_results_objects(
        task, flowcell, target_set
    )
    # Get the barcodes from this dataframe, only one value for each barcode
    barcodes = targets_df["barcode_name"].unique()
    # initialise the NCBITaxa from ete3 for species lookup
    ncbi = NCBITaxa()

    already_created_barcodes = list(
        MappingResult.objects.filter(task=task)
        .values_list("barcode_name", flat=True)
        .distinct()
    )

    # for each barcode in the target_df
    for barcode in barcodes:
        # Check if barcode has mapping result objects created for it already
        if barcode not in already_created_barcodes:
            # for each target, create a Mapping result
            for target in targets:
                # get the tax_id for this target
                tax_id = ncbi.get_name_translator([target])[target][0]
                MappingResult(
                    task=task,
                    barcode_name=barcode,
                    species=target,
                    tax_id=tax_id,
                ).save()

    del ncbi
    # Deep copy the target reads dataframe so changes aren't reflected in the original
    update_num_matches_df = name_df.copy(deep=True)
    # If there's data, update the num_matches and sum_unique fields on mapping results objects for matching rows
    if not update_num_matches_df.empty:
        calculate_num_matches_update(
            update_num_matches_df,
            task,
            num_matches_targets_barcoded_df,
            class_per_bar,
        )

    # Initialise red_df for red alert target reads
    red_reads_df = pd.DataFrame()
    # Get one group for each name
    gb = targets_df.groupby(["name"], as_index=False)
    # for the name of the group and the group itself
    for name, group in gb:
        # append the return value of map_all_the_groups to red reads_df, from each species in group
        red_reads_df = red_reads_df.append(
            map_all_the_groups(
                group,
                name,
                flowcell,
                target_regions_df,
                targets_df,
                task,
                fastas,
            )
        )

    # If none of the reads mapped to inside the groups
    if red_reads_df.empty:
        logger.info(
            "Flowcell id: {} - No reads mapped to dangerous areas!".format(
                flowcell.id
            )
        )
        return

    # results df contains all centrifuge details on reads that mapped to target areas
    results_df = pd.merge(
        targets_df,
        red_reads_df,
        how="inner",
        left_on="read_id",
        right_on="query_seq_name",
    )
    # Deep copy the results df, so changes not reflected in original , one df for barcoded, one df for all reads
    results_df_bar = results_df.copy(deep=True)
    # group non barcoded df by name
    mapping_gb_nb = results_df.groupby("name")
    # set index on non barcoded dataframe to name so we can map results back onto it
    results_df.set_index("name", inplace=True)
    # Get the number of matches for each species as the size of each group
    results_df["red_num_matches"] = mapping_gb_nb.size()
    # Set the barcode name to all reads, so results are irrespective of barcodes
    results_df["barcode_name"] = "All reads"
    # reset index
    results_df.reset_index(inplace=True)
    # Drop duplicates on name, so we only have on value for each species
    results_df.drop_duplicates(subset="name", inplace=True)

    if "name" in results_df.keys():
        results_df.rename(columns={"name": "species"}, inplace=True)

    # For each row, update the red reads value for the mappingTargets entry for any species in the dataframe,
    #  in all reads barcode
    results_df.apply(update_mapped_red_values, args=(task,), axis=1)

    # Barcoded red reads
    # Group by the name and the barcode, so split by each species and the barcode inside that split
    mapping_gb = results_df_bar.groupby(["name", "barcode_name"])
    # Set the index to multi Index on name and barcode_name, so gb calculations map back onto the index
    results_df_bar.set_index(["name", "barcode_name"], inplace=True)

    # TODO this is where we would combine them so the proportion is correct
    # The size of the group is the same as the group size
    results_df_bar["red_num_matches"] = mapping_gb.size()
    # reset the index
    results_df_bar.reset_index(inplace=True)
    # rename the name column to species
    results_df_bar.rename(columns={"name": "species"}, inplace=True)
    # Drop duplicates so a species is unique inside each barcode
    results_df_bar.drop_duplicates(
        subset=["species", "barcode_name"], inplace=True
    )
    # Update the red reads value for the MappingTargets object for the species and barcodes inside this iteration
    results_df_bar.apply(update_mapped_red_values, args=(task,), axis=1)


def create_donut_data_models(row, task):
    """
    Create and return a series of model objects for bulk creation
    :param row: The row of the dataframe
    :type row: pandas.core.series.Series
    :param task: The task objects for this analysis
    :return: A model objects to save into the dataframe
    """

    return DonutData(
        task=task,
        num_matches=row["num_matches"],
        sum_unique=row["sum_unique"],
        name=row["name"],
        tax_rank=row["tax_rank"],
        barcode_name=row["barcode_name"],
        latest=row["latest"],
    )


def calculate_donut_data(df, lineages_df, flowcell, task, tax_rank_filter):
    """
    Calculate the donut data for the donut chart and table
    :param df: The dataframe of centrifuge output, with barcodes as well as all reads, post calculations
    :param lineages_df: The dataframe of lineages for centrifuge output species in this analysis iteration
    :param task: The task object for this analysis
    :param flowcell: The flowcell
    :param tax_rank_filter: A list of the taxonomic ranks in order from superkingdom to species
    :return:
    """
    # Merge the lineages dataframe onto the centrifuge output dataframe so each species has it's lineage
    data_df = pd.merge(df, lineages_df, left_on="tax_id", right_index=True)
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
        "Flowcell id: {} - Bulk inserting new species donut data".format(
            flowcell.id
        )
    )

    donut_to_create_df = donut_df

    # Call output parser to insert the results dataframe
    output_parser(task, donut_to_create_df, "donut", None)


def create_centrifuge_models(row, classified_per_barcode, first):
    """
    Append a CentOutput object to a list for each row in the centrifuge output
    :param row: the row from the data frame
    :param classified_per_barcode: The number of reads classified for each barcode as a dictionary
    :param first: Boolean saying if this is the first iteration, if it isn't a pandas merge adds a suffix onto the code
    :return: The list of newly created objects

    """
    # If there is a NaN in the proportion of classified, redo the calculation
    if (
        row["proportion_of_classified"] == "Unclassified"
        or row["proportion_of_classified"] == np.nan
    ):
        row["proportion_of_classified"] = round(
            row["num_matches"] / classified_per_barcode[row["barcode_name"]], 3
        )

    if first:
        return CentrifugeOutput(
            name=row["name"],
            tax_id=row["tax_id"],
            task=row["task"],
            num_matches=row["num_matches"],
            sum_unique=row["sum_unique"],
            barcode_name=row["barcode_name"],
            proportion_of_classified=row["proportion_of_classified"],
            superkingdom=row["superkingdom"],
            phylum=row["phylum"],
            classy=row["class"],
            order=row["order"],
            family=row["family"],
            genus=row["genus"],
            species=row["species"],
            latest=row["latest"],
        )
    else:

        return CentrifugeOutput(
            name=row["name_x"],
            tax_id=row["tax_id"],
            task=row["task"],
            num_matches=row["num_matches"],
            sum_unique=row["sum_unique"],
            barcode_name=row["barcode_name"],
            proportion_of_classified=row["proportion_of_classified"],
            superkingdom=row["superkingdom_x"],
            phylum=row["phylum_x"],
            classy=row["classy_x"],
            order=row["order_x"],
            family=row["family_x"],
            genus=row["genus_x"],
            species=row["species_x"],
            latest=row["latest"],
        )


def output_parser(task, new_data_df, donut_or_output, metadata):
    """
    Parse the centrifuge output or Donut data, update the values and save into a site readable format
    :param task: The reference to the JobMaster for this metagenomics task
    :param new_data_df: The donut or CentrifugeOutput dataframe calculated this iteration
    :param donut_or_output: Whether the data type is DonutData or CentrifugeOutput
    :param metadata:
    :return:
    """

    metagenomics_task = task
    # If CentrifugeOutput
    if donut_or_output is "output":
        # Get the results from the last Iteration, using the iteration_count (Iteration count) value from the jobmaster
        parsed_data = CentrifugeOutput.objects.filter(
            task=metagenomics_task, latest=task.iteration_count
        )
        # This is not for a donut
        donut = False
        # If DonutData
    else:
        # last Iteration Donut data values
        parsed_data = DonutData.objects.filter(
            task=metagenomics_task, latest=task.iteration_count
        )
        # This is for a donut!
        donut = True
    # Create a dataframe of the previous values
    parsed_data_df = pd.DataFrame(list(parsed_data.values()))
    # class throws a key error cause of poor choices I made, so rename it classy
    if "class" in new_data_df.keys():
        new_data_df.rename(columns={"class": "classy"}, inplace=True)
    # If we have previous data
    if not parsed_data_df.empty:
        # Add one to the iteration count
        new_iteration_count = task.iteration_count + 1
        # If this is for CentrifugeOutput
        if not donut:
            # Merge the lsat Iterations values with this Iterations values
            merged_data_df = pd.merge(
                parsed_data_df,
                new_data_df,
                on=["barcode_name", "tax_id"],
                how="outer",
            )

            merged_data_df["superkingdom_x"].fillna(
                merged_data_df["superkingdom_y"], inplace=True
            )
            merged_data_df["phylum_x"].fillna(
                merged_data_df["phylum_y"], inplace=True
            )
            merged_data_df["classy_x"].fillna(
                merged_data_df["classy_y"], inplace=True
            )
            merged_data_df["order_x"].fillna(
                merged_data_df["order_y"], inplace=True
            )
            merged_data_df["family_x"].fillna(
                merged_data_df["family_y"], inplace=True
            )
            merged_data_df["genus_x"].fillna(
                merged_data_df["genus_y"], inplace=True
            )
            merged_data_df["species_x"].fillna(
                merged_data_df["species_y"], inplace=True
            )
            merged_data_df["name_x"].fillna(
                merged_data_df["name_y"], inplace=True
            )
        # If this is for a donut
        else:
            # Merge the lsat Iterations values with this Iterations values
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
        # Fill those NaNs with ) so we can combine the columns
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
        new_iteration_count = 1
        # Group by barcode name
        gb_bc = new_data_df.groupby("barcode_name")
        # The number of matches per barcode as a dictionary
        classed_per_barcode = gb_bc["num_matches"].sum().to_dict()

        to_save_df = new_data_df

    def divd(row, cl_bar):
        """
        Calculate_proportion_of_classified
        :param row: The df row
        :param cl_bar: The number of reads classified in a barcode, in dict form keyed to the barcode_name
        :return:
        """
        return round(
            (row["num_matches"] / cl_bar[row["barcode_name"]]) * 100, 4
        )

    # Set latest to the updated iteration count
    to_save_df["latest"] = new_iteration_count

    to_save_df["task"] = metagenomics_task

    # If this is for CentrifugeOutput
    if not donut:
        to_save_df["proportion_of_classified"] = to_save_df.apply(
            divd, args=(classed_per_barcode,), axis=1
        )

        to_save_df.rename(columns={"classy": "class"}, inplace=True)

        to_save_df["proportion_of_classified"].fillna(
            "Unclassified", inplace=True
        )

        to_save_df.fillna("Unclassified", inplace=True)

        # The number of reads we have any form of classification for
        reads_classified = to_save_df[
            to_save_df["tax_id"].ne(0)
            & to_save_df["barcode_name"].eq("All reads")
        ]["num_matches"].sum()
        # The number of reads we have completely failed to classify
        reads_unclassified = to_save_df[
            to_save_df["tax_id"].eq(0)
            & to_save_df["barcode_name"].eq("All reads")
        ]["num_matches"].sum()
        # save the values
        metadata.classified = reads_classified
        metadata.unclassified = reads_unclassified
        metadata.save()

        if parsed_data_df.empty:
            # Create a series of CentrifugeOutput objects, one for each row in the dataframe
            to_bulk_save_series = to_save_df.apply(
                create_centrifuge_models,
                args=(classed_per_barcode, True),
                axis=1,
            )
        else:
            to_bulk_save_series = to_save_df.apply(
                create_centrifuge_models,
                args=(classed_per_barcode, False),
                axis=1,
            )
        # Bulk create all those objects in the series
        CentrifugeOutput.objects.bulk_create(
            to_bulk_save_series.values.tolist(), batch_size=1000
        )
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
    # Return the iteration count
    return new_iteration_count


def run_centrifuge(flowcell_job_id):
    """

    Returns nothing.
    -------
     Run the centrifuge command, and open the subprocess to run it. Keep the process open and pipe the reads into
     it Write the reads incrementally to the processes stdin, capture the stdout and parse it into a useful format.
     Once parsed it's deposited into the database, updating any existing records. The organising key is a uuid,
     stored in self.foreign key.
    """

    # The JobMaster object
    task = JobMaster.objects.get(pk=flowcell_job_id)
    # The number of reads for this iteration
    chunk_size = 10000
    # The flowcell the reads are from
    flowcell = task.flowcell

    logger.info(
        "Flowcell id: {} - Running centrifuge on flowcell {}".format(
            flowcell.id, flowcell.name
        )
    )
    logger.info(
        "Flowcell id: {} - job_master_id {}".format(flowcell.id, task.id)
    )
    logger.info(
        "Flowcell id: {} - last_read {}".format(flowcell.id, task.last_read)
    )
    logger.info(
        "Flowcell id: {} - read_count {}".format(flowcell.id, task.read_count)
    )

    # The path to the centrifuge executable
    centrifuge_path = get_env_variable("MT_CENTRIFUGE")
    # The path to the Centrifuge Index
    index_path = get_env_variable("MT_CENTRIFUGE_INDEX")
    # The command to run centrifuge
    cmd = "perl " + centrifuge_path + " -f --mm -k 1 -x " + index_path + " -"

    logger.info("Flowcell id: {} - {}".format(flowcell.id, cmd))

    # Counter for total centOut
    total_centrifuge_output = 0

    # Instance of the Ncbi taxa class, for taxonomic id manipulation
    ncbi = NCBITaxa()

    runs = flowcell.runs.all()

    fastq_df_barcode, last_read, read_count, fasta_objects = call_fetch_reads_cent(
        runs, chunk_size, task.last_read
    )

    if read_count == 0:
        # task.complete = True
        logger.info(
            "Flowcell id: {} - Found 0 reads in the database, for this flowcell. Aborting...".format(
                flowcell.id, chunk_size
            )
        )
        task = JobMaster.objects.filter(pk=flowcell_job_id).first()
        if task.first() is None:
            logger.error("Task appears to have been deleted.")
            pass
        task.running = False
        task.save()
        return

    logger.info(
        "Flowcell id: {} - number of reads found {}".format(
            flowcell.id, read_count
        )
    )

    try:
        # Get the metadata object. Contains the start time, end time and runtime of the task
        metadata = Metadata.objects.get(task=task)

    except (Metadata.DoesNotExist, Metadata.MultipleObjectsReturned):

        metadata = Metadata(task=task)

        metadata.save()

    # Create a fastq string to pass to Centrifuge

    fastq_df_barcode["fasta"] = (
        ">read_id="
        + fastq_df_barcode["read_id"]
        + ",barcode="
        + fastq_df_barcode["barcode_name"]
        + "\n"
        + fastq_df_barcode["sequence"]
    )

    fastqs_data = "\n".join(list(fastq_df_barcode["fasta"]))
    # fastqs_data = "".join([str(">read_id=" + c[0] + ",barcode=" + c[1] + "\n" + c[2] + "\n") for c in fastqs_list
    #                       if c[2] is not None])

    logger.info(
        "Flowcell id: {} - Loading index and Centrifuging".format(flowcell.id)
    )

    # Write the generated fastq file to stdin, passing it to the command
    # Use Popen to run the centrifuge command
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
        logger.info(
            "Flowcell id: {} - Centrifuge error!! {}".format(flowcell.id, err)
        )

    # out is a bytestring so it needs decoding
    centrifuge_output = out.decode()

    # total number of lines of centrifuge output dealt with
    total_centrifuge_output += centrifuge_output.count("\n") - 1

    logger.info(
        "Flowcell id: {} - number of centrifuge output lines is {}".format(
            flowcell.id, total_centrifuge_output
        )
    )

    # output fields is the column headers for the pandas data frame
    output_fields = ["readID", "seqID", "taxID", "numMatches"]

    # create the DataFrame from the output
    df = pd.read_csv(
        StringIO(centrifuge_output), sep="\t", usecols=output_fields
    )
    # If there is no output the dataframe will be empty
    if df.empty:
        logger.info(
            "Flowcell id: {} - No reads found or no centrifuge output."
            " Check above for error".format(flowcell.id)
        )

        task.running = False
        task.save()
        return

    # create DataFrame column unique, where if the number of matches is 1, unique is 1, if matches is more than
    # 1 unique is 0, using numpy condition
    df["unique"] = np.where(
        df["numMatches"] == 1, 1, 0  # condition  # True  # False
    )
    # split out the barcode_name from the readID string
    df = pd.concat(
        [df, df["readID"].str.split(",").str[1].str.split("=").str[1]], axis=1
    )
    # Add the extra column name
    df.columns = [
        "readID",
        "seqID",
        "tax_id",
        "numMatches",
        "unique",
        "barcode_name",
    ]
    # split out the read_id from the readID string
    df = pd.concat(
        [df, df["readID"].str.split(",").str[0].str.split("=").str[1]], axis=1
    )
    # Add the column name
    df.columns = [
        "readID",
        "seqID",
        "tax_id",
        "numMatches",
        "unique",
        "barcode_name",
        "read_id",
    ]

    # get np array of the taxids
    taxid_list = np.unique(df["tax_id"].values)

    # use the ncbi thing to get the species names for the tax ids in the tax ids list
    taxid2name_df = pd.DataFrame.from_dict(
        ncbi.get_taxid_translator(taxid_list), orient="index", columns=["name"]
    )

    # insert the taxid2name dict items into the dataframe name columns
    df = pd.merge(
        df, taxid2name_df, how="outer", left_on="tax_id", right_index=True
    )

    # any NaNs replace with Unclassified
    df["name"].fillna("Unclassified", inplace=True)

    #  set the index to be taxID
    df = df.set_index("tax_id")
    # group by the barcode
    gb_bc = df.groupby(["barcode_name"])

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
    barcodes = list(df["barcode_name"].unique())
    # Add an all reads string to the list
    barcodes.append("All reads")
    # Get the reads classified in each barcode into a lookup dictionary
    #  so we can work out proportions inside the barcode
    classified_per_barcode = {
        barcode: CentrifugeOutput.objects.filter(
            task=task, barcode_name=barcode
        ).aggregate(Sum("num_matches"))["num_matches__sum"]
        for barcode in barcodes
    }
    # If this is the first iteration, classified per barcode will be NaN
    if classified_per_barcode["All reads"] is None:
        # Set it to 0
        classified_per_barcode["All reads"] = 0
    # for the name and group_df in the barcode group by object
    for name, group in gb_bc:
        # If there is no value for that barcode
        if classified_per_barcode[name] is None:
            classified_per_barcode[name] = 0
        # set the value
        classified_per_barcode[name] += group.shape[0]

    # Unique targets in this set of regions
    temp_targets = (
        MappingTarget.objects.filter(target_set=task.target_set)
        .values_list("species", flat=True)
        .distinct()
    )

    logger.info(
        "Flowcell id: {} - Targets are {}".format(flowcell.id, temp_targets)
    )
    # Get the target reads rows in a seperate dataframe
    name_df = df[df["name"].isin(temp_targets)]

    logger.info(
        "Flowcell id: {} - The dataframe target dataframe is {}".format(
            flowcell.id, name_df.head(n=20)
        )
    )
    # initialise the barcode_df
    barcode_df = pd.DataFrame()

    # TODO vectorise these bad boys
    logger.info(
        "Flowcell id: {} - The dataframe shape is {}".format(
            flowcell.id, df.shape
        )
    )

    df2 = df.reset_index()
    # Drop duplicates in the dataframe
    df2 = df2.drop_duplicates(keep="first")
    # Make this the results for all reads in the centrifuge output
    df2["barcode_name"] = "All reads"

    # Include these reads in the lookup dict
    classified_per_barcode["All reads"] += df2.shape[0]

    del df2
    # Create a dataframe with the calculations applied to each barcode
    barcode_df = gb_bc.apply(
        calculate_barcoded_values, barcode_df, classified_per_barcode
    )
    # Reset the index
    barcode_df.reset_index(inplace=True)
    # rename the columns to useful names
    barcode_df.rename(
        columns={"unique": "sum_unique", 0: "num_matches"}, inplace=True
    )
    # Set the tax_id as the index
    barcode_df.set_index("tax_id", inplace=True)
    # TODO this may be slow
    # Get the previous Centrifuge Output objects so we can get there num matches and sum unique
    df_data = CentrifugeOutput.objects.filter(
        species__in=temp_targets, task=task, barcode_name__in=barcodes
    ).values("species", "tax_id", "num_matches", "barcode_name", "sum_unique")
    # Get the results into  a dataframe
    num_matches_per_target_df = pd.DataFrame(list(df_data))

    logger.info(
        "Flowcell id: {} - The previous number of matches dataframe is {}".format(
            flowcell.id, num_matches_per_target_df
        )
    )

    already_created_barcodes = list(
        MappingResult.objects.filter(task=task)
        .values_list("barcode_name", flat=True)
        .distinct()
    )

    # for each barcode in the target_df
    for barcode in barcodes:
        # Check if barcode has mapping result objects created for it already
        if barcode not in already_created_barcodes:
            # for each target, create a Mapping result
            for target in temp_targets:
                # get the tax_id for this target
                tax_id = ncbi.get_name_translator([target])[target][0]
                mr = MappingResult(
                    task=task,
                    barcode_name=barcode,
                    species=target,
                    tax_id=tax_id,
                )
                mr.save()
    # Call map the reads on the targets dataframe

    map_the_reads(
        name_df,
        task,
        flowcell,
        num_matches_per_target_df,
        temp_targets,
        classified_per_barcode,
        fasta_objects,
    )

    # delete these columns, no longer needed
    df.drop(
        columns=[
            "readID",
            "seqID",
            "numMatches",
            "unique",
            "barcode_name",
            "read_id",
        ],
        inplace=True,
    )

    # remove duplicate rows in the data frame,so we only have one entry for each species
    df.drop_duplicates(keep="first", inplace=True)

    # filter out no ranks, like root and cellular organism
    tax_rank_filter = [
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]

    # drop all duplicated lines to keep only one for each entry, creating atomicity
    df.reset_index(inplace=True)

    df = df.drop_duplicates(keep="first")

    # Make this the results for all reads in the centrifuge output
    df["barcode_name"] = "All reads"

    df.set_index("tax_id", inplace=True)
    # Add the barcode dataframe onto the dataframe, so now we have all reads and barcodes
    df = df.append(barcode_df, sort=True)
    # Set the flowcell object into a series, broadcast to each row
    df["flowcell"] = flowcell

    # same with task from Jobmaster
    # Get the task again so we don't over write the paused status
    task = JobMaster.objects.get(pk=task.id)
    df["task"] = task

    df.reset_index(inplace=True)
    # Create the lineages df, a df with the lineage of each species in the full results dataframe
    lineages_df = insert_new_lineages(ncbi, df, tax_rank_filter, flowcell)

    # Calculate the donut data, and it is inserted in this function
    try:
        calculate_donut_data(df, lineages_df, flowcell, task, tax_rank_filter)
    except KeyError as e:
        logger.info("PROBLEMS {}".format(e))
    #        return

    cent_to_create_df = df

    logger.info(
        "Flowcell id: {} - Centrifuge to lineages dataframe {}".format(
            flowcell.id, lineages_df.head()
        )
    )
    # Merge the lineages dataframe onto the centrifuge dataframe, so each row has the lineage of the species of that row
    cent_to_create_df = pd.merge(
        cent_to_create_df,
        lineages_df,
        how="inner",
        left_on="tax_id",
        right_index=True,
    )
    # Broadcast the taskobject down the dataframe
    cent_to_create_df["task"] = task
    # Parse that dataframe, and store the results, returning the new latest iteration count
    new_latest = output_parser(task, cent_to_create_df, "output", metadata)
    # Get the task again so we don't over write the paused status
    task = JobMaster.objects.get(pk=task.id)
    # Set the iteration count on the job master
    task.iteration_count = new_latest
    # Save the new value
    task.save()
    # Delete the values from the old iteration for the Centrifuge Output data
    CentrifugeOutput.objects.filter(task=task, latest__lt=new_latest).delete()
    # Delete the values from the old iteration from the Donut Data
    DonutData.objects.filter(task=task, latest__lt=new_latest).delete()

    # Update the metadata for the headers

    # Get the start time
    start_time = metadata.start_time.replace(tzinfo=None)
    # Get the end time
    end_time = timezone.now().replace(tzinfo=None)
    # Calculate the finish time
    metadata.finish_time = str(end_time - start_time)
    # Save the values
    metadata.save()
    # Refresh the task
    task = JobMaster.objects.get(pk=task.id)
    # Set running to false
    task.running = False
    # OK - this returns the last record of all fastqs, not the last record of the chunk.
    task.last_read = last_read
    # Set the read count
    task.read_count = task.read_count + read_count
    # Save the new task values
    task.save()

    logger.info(
        "Flowcell id: {} - New last_read_id is - {}".format(
            flowcell.id, task.last_read
        )
    )
    logger.info(
        "Flowcell id: {} - Total CentOut lines - {}".format(
            flowcell.id, total_centrifuge_output
        )
    )
    logger.info("Flowcell id: {} - Finished!".format(flowcell.id))
