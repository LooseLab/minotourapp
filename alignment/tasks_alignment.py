from __future__ import absolute_import, unicode_literals

import os
import subprocess
from collections import defaultdict
from io import StringIO

import numpy as np
import pandas as pd
from celery import task
from celery.utils.log import get_task_logger
from django.conf import settings
from django.db.models import Max

from alignment.models import PafRoughCov, PafSummaryCov
from reads.models import FastqRead, JobMaster
from reference.models import ReferenceInfo

# Set up the logger to write to logging file
logger = get_task_logger(__name__)


@task()
def run_minimap2_alignment(job_master_id, streamed_reads=[]):
    """
    Run a minimap2 based alignment.
    Parameters
    ----------
    job_master_id: int
        The Pk of the JobMaster entry in the database
    streamed_reads: list of dict
        A list of read dictionaries, if reads are being streamed straight into the task
    Returns
    -------
    None

    """
    # Get the JobMaster
    job_master = JobMaster.objects.get(pk=job_master_id)
    if not job_master.reference:
        raise ValueError("Missing Reference file. Please sort out.")
    flowcell = job_master.flowcell
    avg_read_length = flowcell.average_read_length
    if avg_read_length == 0:
        logger.error(
            f"Average read length is zero Defaulting to 450, but this is an error."
        )
        avg_read_length = 450
    if not streamed_reads:
        # The chunk size of reads from this flowcell to fetch from the database
        # aim for 50 megabases
        desired_yield = 50 * 1000000
        chunk_size = round(desired_yield / avg_read_length)
        logger.debug(f"Fetching reads in chunks of {chunk_size} for alignment.")
        fasta_objects = FastqRead.objects.filter(
            flowcell_id=flowcell.id, id__gt=job_master.last_read
        )[:chunk_size]
        fasta_df_barcode = pd.DataFrame().from_records(
            fasta_objects.values(
                "read_id",
                "barcode_name",
                "sequence",
                "id",
                "run_id",
                "type__name",
                "type_id",
                "run_id",
                "barcode_id",
                "is_pass",
            )
        )
        if fasta_objects:
            last_read = fasta_df_barcode.tail(1).iloc[0].id
        if fasta_df_barcode.shape[0] == 0:
            logger.debug("No fastq found. Skipping this iteration.")
            job_master.running = False
            job_master.save()
            return
    else:
        last_read = job_master.last_read
        fasta_df_barcode = pd.DataFrame(streamed_reads)
        fasta_df_barcode = fasta_df_barcode[fasta_df_barcode["is_pass"]]
        fasta_objects = streamed_reads

    read_count = fasta_df_barcode.shape[0]
    # The chunk size of reads from this flowcell to fetch from the database

    logger.debug(
        f"Flowcell id: {flowcell.id} - Running minimap2 on flowcell {flowcell.name}"
    )
    logger.debug(f"Flowcell id: {flowcell.id} - job_master_id {job_master.id}")
    logger.debug(f"Flowcell id: {flowcell.id} - reference {job_master.reference.name}")
    logger.debug(f"Flowcell id: {flowcell.id} - last_read {job_master.last_read}")
    logger.debug(f"Flowcell id: {flowcell.id} - read_count {job_master.read_count}")
    logger.debug(f"Flowcell id: {flowcell.id} - number of reads found {read_count}")
    # If we have pulled back reads, call fasta
    if read_count > 0:
        align_reads(job_master.id, fasta_df_barcode)

    # Update the JobMaster with the metadata after finishing this iteration
    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = False
    job_master.last_read = last_read
    job_master.read_count = job_master.read_count + read_count
    job_master.save()


def save_paf_store_summary(row, job_master_id):
    """
    Save the paf store summaries
    Parameters
    ----------
    row: pd.core.series.Series
        The row of the dataframe we are saving
    job_master_id: int
        The primary key of the databases task record

    Returns
    -------

    """

    paf_summary_cov, created = PafSummaryCov.objects.get_or_create(
        job_master_id=job_master_id,
        barcode_name=row["barcode_name"],
        barcode_id=row["barcode_id"],
        chromosome_id=row["chromosome_pk"],
        reference_line_length=row[6],
        read_type_id=row["type_id"],
    )

    paf_summary_cov.total_yield += row["yield"]
    paf_summary_cov.read_count += row["read_count"]
    paf_summary_cov.coverage = round(paf_summary_cov.total_yield / row[6])
    paf_summary_cov.average_read_length = round(
        paf_summary_cov.total_yield / paf_summary_cov.read_count, 2
    )
    paf_summary_cov.save()


def generate_paf_summary_cov_objects(row, flowcell_pk, job_master_pk, reference_pk):
    """
    Generate Paf summary coverage ORM objects on a dataframe row
    Parameters
    ----------
    row: pd.core.series.Series
        The row of a dataframe
    flowcell_pk: int
        The primary key of the flowcell database record
    job_master_pk: int
        The primary key of the job master database record
    reference_pk: int
        The reference primary key database record


    Returns
    -------
    alignment.models.PafSummaryCov

    """
    return PafRoughCov(
        job_master_id=job_master_pk,
        run_id=row["run_id"],
        flowcell_id=flowcell_pk,
        read_type_id=row["type_id"],
        is_pass=row["is_pass"],
        barcode_id=row["barcode_id"],
        reference_id=reference_pk,
        chromosome_id=row["chromosome_pk"],
        bin_position_start=row["left_bin"],
        bin_position_end=row["right_bin"],
        bin_coverage=row["bin_coverage"],
    )


def paf_rough_coverage_calculations(df, job_master, max_chromosome_length):
    """
    Calculate the rough coverage binned by 10 base intervals in the paf, then dump them into database
    Parameters
    ----------
    df: pandas.core.frame.DataFrame
        The pandas dataframe containing our paf output and the barcode, read type and chromosome name
    job_master: reads.models.JobMaster
        The primary key of the job master
    max_chromosome_length: int
        The maximum length out of reference contig chromosome

    Returns
    -------

    """
    bin_width = 10
    bins = [i for i in range(0, max_chromosome_length + bin_width, bin_width)]
    df["cut"] = pd.cut(df[7], bins)
    df["left_bin"] = df["cut"].apply(lambda x: x.left)
    df["right_bin"] = df["cut"].apply(lambda x: x.right)
    if "barcode_name" not in df.columns:
        df.reset_index(inplace=True)
    df.set_index(["type__name", "tsn", "barcode_name", "cut"], inplace=True)
    df["bin_coverage"] = df.groupby(["type__name", "tsn", "barcode_name", "cut"]).size()
    df = df.loc[~df.index.duplicated()]
    df["model_objects"] = df.apply(
        generate_paf_summary_cov_objects,
        axis=1,
        args=(job_master.flowcell.id, job_master.id, job_master.reference.id),
    )

    PafRoughCov.objects.bulk_create(df["model_objects"].values.tolist(), batch_size=5000)


def paf_summary_calculations(df, job_master_id):
    """
    Do dataframe aggregations on a per chromosome, barcode and readtype combination level
    and save the resulting pafstore
    Parameters
    ----------
    df: pd.core.frame.DataFrame
        Merged dataframe of paf output and fastq dataframe
    job_master_id: int
        The primary key of this task database record

    Returns
    -------

    """
    # TODO read_id overwritten by index change before function call somewhere
    if "type__name" in df.keys():
        df.set_index(["type__name", "tsn", "barcode_name"], inplace=True)
    df[["yield", "read_count"]] = df.groupby(
        level=["type__name", "tsn", "barcode_name"]
    ).agg({"mapping_length": np.sum, 1: np.size})
    df = df.loc[~df.index.duplicated()]
    df.reset_index(inplace=True)
    df.apply(save_paf_store_summary, args=(job_master_id,), axis=1)


def align_reads(job_master_id, fasta_df_barcode):
    """
    Call minimap2 to align reads to the reference sequence listed on the JobMaster.

    Parameters
    ----------
    fastas: list
        A list of FastqRead objects.
    job_master_id: int
        The PK of the Jobmaster entry.
    fasta_df_barcode: pandas.core.frame.DataFrame
        DataFrame containing all the sequence data that we need.

    Returns
    -------
    last_read: int
        The ID of the last read that we dealt with.
    """
    # The on disk directory where we store the references
    reference_location = getattr(settings, "REFERENCE_LOCATION", None)
    # The location of the mimimap2 executable
    minimap2 = getattr(settings, "MINIMAP2", None)
    # The JobMaster object
    job_master = JobMaster.objects.get(pk=job_master_id)
    # Unpack metadata about this job
    last_read = job_master.last_read
    reference_pk = job_master.reference.id
    flowcell = job_master.flowcell
    #
    if not minimap2:
        logger.error("Can not find minimap2 executable - stopping task.")
        return
    # We're starting fresh boys
    if last_read is None:
        last_read = 0
    # The debugrmation about the reference for this alignment
    reference_info = ReferenceInfo.objects.get(pk=reference_pk)
    # A dictionary with all the results for each chromosome
    chromosome_dict = {}
    # Get all the chromosomes for this reference
    chromosomes = reference_info.referencelines.all()
    max_chromosome_length = chromosomes.aggregate(max_length=Max("chromosome_length"))[
        "max_length"
    ]
    # For each Chromosome in the chromosomes, set the name of the chromosome as key, and the object as the value
    for chromosome in chromosomes:
        chromosome_dict[chromosome.line_name] = chromosome.id
    # Create the series that contains the read_id and sequence as a correctly formatted fasta string
    fasta_df_barcode["fasta"] = (
        ">" + fasta_df_barcode["read_id"] + "\n" + fasta_df_barcode["sequence"]
    )
    # Create one string formatted correctly for fasta with input data
    fasta_data = "\n".join(list(fasta_df_barcode["fasta"]))
    # Create the minimap command
    # If we have a minimap2 index we want to use it.
    if reference_info.minimap2_index_file_location:
        reference_file = reference_info.minimap2_index_file_location
    else:
        reference_file = reference_info.filename
    cmd = f"{minimap2} -x map-ont -t 4 --secondary=no {os.path.join(reference_location, reference_file)} -"
    logger.debug(f"Flowcell id: {flowcell.id} - Calling minimap - {cmd}")
    # Use subprocess to run minimap
    proc = subprocess.Popen(
        cmd.split(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdin=subprocess.PIPE,
        universal_newlines=True,
    )
    # Run the command and get a tuple of output
    (out, err) = proc.communicate(input=fasta_data)
    if err:
        logger.error(err)
    logger.debug(f"Flowcell id: {flowcell.id} - Finished minimap - {cmd}")
    # Get the results of the Paf output
    paf = out
    paf_line_count = paf.count("n")
    logger.debug(f"Flowcell id: {flowcell.id} - Found {paf_line_count} paf records")
    # Store the data to become a pafsumamrycov record in DB
    paf_summary_cache = defaultdict(dict)
    paf_df = pd.read_csv(StringIO(out), sep="\t", header=None)
    paf_df = paf_df.rename(columns={5: "tsn"})
    paf_df["chromosome_pk"] = paf_df["tsn"].map(chromosome_dict)
    paf_df.set_index(0, inplace=True)
    paf_df = paf_df[paf_df[11] >= 11]
    fasta_df_barcode.set_index("read_id", inplace=True)
    df = pd.merge(
        paf_df,
        fasta_df_barcode[
            ["barcode_name", "type__name", "type_id", "barcode_id", "is_pass", "run_id"]
        ],
        how="left",
        left_index=True,
        right_index=True,
    )
    df["mapping_length"] = df[8] - df[7]
    if "no barcode" in fasta_df_barcode["barcode_name"].unique():
        df["barcode_name"] = "All reads"
        paf_summary_calculations(df, job_master_id)
    else:
        paf_summary_calculations(df, job_master_id)
        df.reset_index(inplace=True)
        df["barcode_name"] = "All reads"
        print("seconds call")
        print(df.head())
        print(df.keys())
        # save results for all reads
        paf_summary_calculations(df, job_master_id)
    paf_rough_coverage_calculations(df, job_master, max_chromosome_length)

    return last_read
