from __future__ import absolute_import, unicode_literals

from collections import namedtuple

import numpy as np
import pandas as pd
import redis
from celery.utils.log import get_task_logger
from django.conf import settings
from django.db.models import F

from alignment.mapper import MAP
from alignment.models import PafRoughCov, PafSummaryCov
# from celery import task, Task
from minotourapp.celery import app
from reads.models import FastqRead, JobMaster
from reference.models import ReferenceInfo

# Set up the logger to write to logging file
logger = get_task_logger(__name__)
redis_instance = redis.StrictRedis(
    host="127.0.0.1", port=6379, db=0, decode_responses=True
)


# TODO don't stream reads, pull them out of redis here pass UUid, add ignore results
# @task(base=MappingServer)
@app.task
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
    avg_read_length = int(flowcell.average_read_length)
    if avg_read_length == 0:
        logger.error(
            f"Average read length is zero Defaulting to 450, but this is an error."
        )
        avg_read_length = 450
    if not streamed_reads:
        # The chunk size of reads from this flowcell to fetch from the database
        # aim for 50 megabases
        desired_yield = 100 * 1000000
        chunk_size = round(desired_yield / avg_read_length)
        logger.debug(f"Fetching reads in chunks of {chunk_size} for alignment.")
        fasta_df_barcode = pd.DataFrame().from_records(
            FastqRead.objects.filter(
                flowcell_id=flowcell.id, id__gt=job_master.last_read
            )[:chunk_size].values(
                "read_id",
                "barcode_name",
                "sequence",
                "id",
                "run_id",
                "type__name",
                "run_id",
                "barcode_id",
                "is_pass",
                read_type_id=F("type_id"),
            )
        )
        if not fasta_df_barcode.empty:
            last_read = fasta_df_barcode.tail(1).iloc[0].id
        if fasta_df_barcode.empty:
            logger.debug("No fastq found. Skipping this iteration.")
            job_master.running = False
            job_master.save()
            return
    else:
        last_read = job_master.last_read
        ## So these two dataframes do not have the same column names.
        fasta_df_barcode = pd.DataFrame(streamed_reads)
        # Todo: this seems to only process pass reads?
        # fasta_df_barcode = fasta_df_barcode[fasta_df_barcode["is_pass"]]
        fasta_df_barcode = fasta_df_barcode.rename(
            columns={"type": "read_type_id", "barcode": "barcode_id"}
        )
        fasta_df_barcode["type__name"] = fasta_df_barcode["read_type_id"]

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
        align_reads_factory(job_master.id, fasta_df_barcode, MAP)

    # Update the JobMaster with the metadata after finishing this iteration
    redis_instance.decr("minimaptasks")
    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = False
    job_master.last_read = last_read
    job_master.read_count = job_master.read_count + read_count
    job_master.save()


def save_summary_coverage(row, job_master_id):
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
        read_type_id=row["read_type_id"],
    )
    paf_summary_cov.total_yield += row["yield"]
    paf_summary_cov.read_count += row["read_count"]
    paf_summary_cov.coverage = round(paf_summary_cov.total_yield / int(row[6]), 2)
    paf_summary_cov.average_read_length = round(
        paf_summary_cov.total_yield / paf_summary_cov.read_count, 2
    )
    paf_summary_cov.save()


def fetch_rough_cov_to_update(job_master_id, df):
    """
    Returns a DataFrame of PafRoughCov bins we have data for
    Parameters
    ----------
    job_master_id: int
        The primary key of the JobMaster record for this task
    df: pandas.core.frame.DataFrame
        List of bin position_starts
    Returns
    -------
    pd.core.frame.DataFrame
    """

    querysets = (
        df.reset_index().groupby("chromosome_pk")[["chromosome_pk", "bin_position_start"]]
        .apply(
            lambda x: PafRoughCov.objects.filter(
                chromosome_id=np.unique(x.chromosome_pk.values)
            ).filter(bin_position_start__in=np.unique(x.bin_position_start))
        )
        .values
    )
    s = pd.Series(querysets).explode().to_frame()
    s.dropna(inplace=True)
    if not s.empty:
        lookup_dict = {rough_cov.id: rough_cov for rough_cov in s[0].values.tolist()}
        df_old_coverage = pd.DataFrame(
            s[0].apply(lambda x: x.__dict__).values.tolist()
        )
    else:
        df_old_coverage = pd.DataFrame()
        lookup_dict = {}
    if not df_old_coverage.empty:
        df_old_coverage.drop(
            columns=["_state", "grouprun_id", "bin_position_end", "is_pass"], inplace=True
        )
        df_old_coverage.rename(columns={"chromosome_id": "chromosome_pk"}, inplace=True)
        df_old_coverage["ORM"] = df_old_coverage["id"].map(lookup_dict)
        df_old_coverage = df_old_coverage.set_index(
            [
                "read_type_id",
                "chromosome_pk",
                "barcode_id",
                "bin_position_start",
            ]
        )
        df_old_coverage = df_old_coverage.loc[~df_old_coverage.index.duplicated()]
    return df_old_coverage


def update_orm_objects(row):
    """
    Update the orm objects for bulk updating

    Parameters
    ----------
    row: pandas.core.series.Series

    Returns
    -------

    """
    row["ORM"].bin_coverage = row["bin_coverage"]
    return row["ORM"]


def generate_paf_rough_cov_objects(row, flowcell_pk, job_master_pk, reference_pk):
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
        read_type_id=row["read_type_id"],
        is_pass=row["is_pass"],
        barcode_id=row["barcode_id"],
        reference_id=reference_pk,
        chromosome_id=row["chromosome_pk"],
        bin_position_start=row["bin_position_start"],
        bin_coverage=row["bin_coverage"],
    )


def paf_rough_coverage_calculations(df, job_master):
    """
    Calculate the rough coverage binned by 10 base intervals in the paf, then dump them into database
    Parameters
    ----------
    df: pandas.core.frame.DataFrame
        The pandas dataframe containing our paf output and the barcode, read type and chromosome name
    job_master: reads.models.JobMaster
        The primary key of the job master

    Returns
    -------

    """
    flowcell_id = int(job_master.flowcell.id)
    reference_id = int(job_master.reference.id)
    job_master_id = int(job_master.id)
    df["mapping_start_bin_start"] = df["mapping_start"].sub(5).div(10).round().mul(10)
    df["mapping_start_bin_end"] = df["mapping_start_bin_start"].add(10)
    df["mapping_end_bin_start"] = df["mapping_end"].sub(5).div(10).round().mul(10)
    df["mapping_end_bin_end"] = df["mapping_end_bin_start"].add(10)
    if "barcode_name" not in df.columns:
        df.reset_index(inplace=True)
    df.set_index(
        ["read_type_id", "chromosome_pk", "barcode_id", "mapping_start_bin_start"],
        inplace=True,
    )
    # df.sort_index(inplace=True)
    df["bin_coverage"] = df.groupby(
        ["read_type_id", "chromosome_pk", "barcode_id", "mapping_start_bin_start"]
    ).size()
    df_mapping_start_bins = df.loc[~df.index.duplicated()][
        ["bin_coverage", "run_id", "is_pass", "mapping_start_bin_end"]
    ]
    df_mapping_start_bins.index = df_mapping_start_bins.index.rename(
        ["read_type_id", "chromosome_pk", "barcode_id", "bin_position_start"]
    )
    df.reset_index(inplace=True)
    df.set_index(
        ["read_type_id", "chromosome_pk", "barcode_id", "mapping_end_bin_start"],
        inplace=True,
    )
    # df.sort_index(inplace=True)
    df["end_bin_coverage"] = df.groupby(
        ["read_type_id", "chromosome_pk", "barcode_id", "mapping_end_bin_start"]
    ).size()
    # end_bins are
    df_mapping_end_bins = df.loc[~df.index.duplicated()][
        ["end_bin_coverage", "run_id", "is_pass", "mapping_start_bin_end"]
    ]
    df_mapping_end_bins = df_mapping_end_bins.rename(
        columns={"end_bin_coverage": "bin_coverage"}
    )
    df_mapping_end_bins["bin_coverage"] = -df_mapping_end_bins["bin_coverage"]
    df_mapping_end_bins.index = df_mapping_end_bins.index.rename(
        ["read_type_id", "chromosome_pk", "barcode_id", "bin_position_start"]
    )
    results = df_mapping_start_bins.append(df_mapping_end_bins).sort_index()
    results["bin_coverage"] = results.groupby(df_mapping_start_bins.index.names[:-1])[
        "bin_coverage"
    ].agg(np.cumsum)
    del df_mapping_end_bins, df_mapping_start_bins
    df_old_coverage = fetch_rough_cov_to_update(job_master_id, results)
    results = results.loc[~results.index.duplicated()]
    if not df_old_coverage.empty:
        df_old_coverage = df_old_coverage.loc[
            df_old_coverage.index.isin(results.index.values)
        ]
        df_old_coverage["bin_coverage"] = (
            df_old_coverage["bin_coverage"]
            + results.loc[
                results.index.isin(df_old_coverage.index.values)
            ]["bin_coverage"]
        )
        df_old_coverage["ORM"] = df_old_coverage.apply(update_orm_objects, axis=1)
        PafRoughCov.objects.bulk_update(
            df_old_coverage["ORM"].values, ["bin_coverage"], batch_size=5000
        )
    results = results.loc[
        ~results.index.isin(df_old_coverage.index.values)
    ]
    if not results.empty:
        results.reset_index(inplace=True)
        results["model_objects"] = results.apply(
            generate_paf_rough_cov_objects,
            axis=1,
            args=(flowcell_id, job_master_id, reference_id),
        )
        a = PafRoughCov.objects.bulk_create(
            results["model_objects"].values.tolist(), batch_size=5000
        )


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
    ).agg({"mapping_length": np.sum, "sequence_length": np.size})
    df[6] = df[6].astype(int)
    df = df.loc[~df.index.duplicated()]
    df.reset_index(inplace=True)
    df.apply(save_summary_coverage, args=(job_master_id,), axis=1)


def align_reads_factory(job_master_id, fasta_df_barcode, super_function):
    """
    Sends reads to an alignment factory to map reads to the reference.
    Parameters
    ----------
    job_master_id: int
        The PK of the Jobmaster entry.
    fasta_df_barcode: pandas.core.frame.DataFrame
        DataFrame containing all the sequence data that we need.
    super_function: MappingServer
        Mapping server class for managing mappy

    Returns
    -------
    last_read_pk: int
        The ID of the last read that we dealt with.
    """
    # The location of the mimimap2 executable
    minimap2 = getattr(settings, "MINIMAP2", None)
    # The JobMaster object
    job_master = JobMaster.objects.get(pk=job_master_id)
    # Unpack metadata about this job
    last_read_pk = int(job_master.last_read)
    reference_pk = int(job_master.reference.id)
    #
    if not minimap2:
        logger.error("Can not find minimap2 executable - stopping task.")
        # return
    # We're starting fresh boys
    if last_read_pk is None:
        last_read_pk = 0
    # The debugrmation about the reference for this alignment
    reference_info = ReferenceInfo.objects.get(pk=reference_pk)
    # A dictionary with all the results for each chromosome
    chromosome_dict = {}
    # Get all the chromosomes for this reference
    chromosomes = reference_info.referencelines.all()
    # For each Chromosome in the chromosomes, set the name of the chromosome as key, and the object as the value
    for chromosome in chromosomes:
        chromosome_dict[chromosome.line_name] = chromosome.id
    ## This gives us a json string that we can query against the alignment server.
    MapSequence = namedtuple("MapSequence", ["read_id", "sequence"])
    query_list = (
        fasta_df_barcode[["read_id", "sequence"]]
        .apply(lambda x: MapSequence(x.read_id, x.sequence), axis=1)
        .values.tolist()
    )

    if not super_function.valid(reference_info.name):
        super_function.add_reference(
            reference_info.name, reference_info.minimap2_index_file_location
        )
    # send the sequence to the mapper
    results = []
    for read_tuple in query_list:
        results.extend(super_function.map_sequence(reference_info.name, read_tuple))

    test_df = pd.DataFrame(results)
    paf_df = test_df[0].str.split("\t", expand=True)
    del test_df
    # return test_df

    paf_df = paf_df.rename(
        columns={0: "read_id", 1: "sequence_length", 5: "tsn", 7: "mapping_start", 8: "mapping_end", 11: "mq"}
    )
    paf_df = paf_df.apply(pd.to_numeric, errors="ignore")

    paf_df["chromosome_pk"] = paf_df["tsn"].map(chromosome_dict)
    paf_df.set_index("read_id", inplace=True)
    paf_df = paf_df[paf_df["mq"] >= 40]
    fasta_df_barcode.set_index("read_id", inplace=True)
    df = pd.merge(
        paf_df,
        fasta_df_barcode[
            [
                "barcode_name",
                "type__name",
                "read_type_id",
                "barcode_id",
                "is_pass",
                "run_id",
            ]
        ],
        how="left",
        left_index=True,
        right_index=True,
    )
    df["mapping_length"] = df["mapping_end"] - df["mapping_start"]
    # return
    # Todo: should this be "No barcode" rather than "no barcode"?
    # TODO should these go to a chain, final in chain deletes the paf results in redis
    if "No barcode" in fasta_df_barcode["barcode_name"].unique():
        df["barcode_name"] = "All reads"
        paf_summary_calculations(df, job_master_id)
        paf_rough_coverage_calculations(df, job_master)
    else:
        paf_summary_calculations(df, job_master_id)
        df.reset_index(inplace=True)
        df2 = df
        paf_rough_coverage_calculations(df2, job_master)
        if "barcode_name" not in df.keys():
            df.reset_index(inplace=True)
        df["barcode_name"] = "All reads"
        # save results for all reads
        paf_summary_calculations(df, job_master_id)
        paf_rough_coverage_calculations(df2, job_master)

    return last_read_pk
