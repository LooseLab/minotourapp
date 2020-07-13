from __future__ import absolute_import, unicode_literals

from collections import namedtuple

import numpy as np
import pandas as pd
import redis
from celery.utils.log import get_task_logger
from django.db.models import F, Max
from numba import njit

from alignment.mapper import MAP
from alignment.models import PafRoughCov, PafSummaryCov
from minotourapp.celery import app
from reads.models import FastqRead, JobMaster
from reference.models import ReferenceInfo

# Set up the logger to write to logging file
logger = get_task_logger(__name__)
redis_instance = redis.StrictRedis(
    host="127.0.0.1", port=6379, db=0, decode_responses=True
)


# TODO don't stream reads, pull them out of redis here pass UUid, add ignore results
@app.task
def remove_old_references():
    """
    Remove old references from the mapping task, if they haven't been used for the pre specified time
    Returns
    -------

    """
    MAP.reference_monitor()


@app.task
def run_minimap2_alignment(job_master_id, streamed_reads=None):
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
        desired_yield = 25 * 1000000
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


def fetch_rough_cov_to_update(job_master_id, positions):
    """
    Returns a DataFrame of PafRoughCov bins we have data for
    Parameters
    ----------
    job_master_id: int
        The primary key of the JobMaster record for this task
    positions: pandas.core.frame.Series
        List of bin position_starts
    Returns
    -------
    pd.core.frame.DataFrame
    """
    logger.debug("\nPositions")
    logger.debug(positions)
    querysets = (
        positions.groupby("chromosome_pk").apply(lambda x: PafRoughCov.objects.filter(
            chromosome_id=x.name, job_master_id=job_master_id, bin_position_start__in=x.values[0]
        )).values
    )
    logger.debug("Querysets:")
    logger.debug(querysets)
    if not querysets[0]:
        logger.debug("No previous bins found")
        df_old_coverage = pd.DataFrame()
        return df_old_coverage
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


@njit
def go_fast2(start_stop_bools, bin_width):
    """
    Expand all possible bins between the given starts and ends for a bin width
    Parameters
    ----------
    start_stop_bools: np.ndarray
    bin_width: np.uint32
        The width of the bins

    Returns
    -------

    """
    print("Now we are creating the positions")
    positions = [np.uint32(0)]
    for start, stop, is_start, next_is_start in start_stop_bools:
        if not is_start and next_is_start:
            continue
        else:
            a = np.arange(start, stop + bin_width, bin_width, dtype=np.uint32)
            for position in a:
                positions.append(position)
            return positions


def paf_rough_coverage_calculations(df, job_master, longest_chromosome):
    """
    Calculate the rough coverage binned by 10 base intervals in the paf, then dump them into database
    Parameters
    ----------
    df: pandas.core.frame.DataFrame
        The pandas dataframe containing our paf output and the barcode, read type and chromosome name
    job_master: reads.models.JobMaster
        The primary key of the job master
    longest_chromosome: int
        The longest chromosome base length

    Returns
    -------

    """
    flowcell_id = int(job_master.flowcell.id)
    reference_id = int(job_master.reference.id)
    job_master_id = int(job_master.id)
    bin_width = 10
    if longest_chromosome / bin_width > 500000:
        bin_width = round(longest_chromosome / 500000 / 10) * 10
    if bin_width < 10:
        bin_width = 10
    rounder = bin_width / 2
    df["mapping_start_bin_start"] = df["mapping_start"].sub(rounder).div(bin_width).round().mul(bin_width)
    df["mapping_end_bin_start"] = df["mapping_end"].sub(rounder).div(bin_width).round().mul(bin_width)
    if "barcode_id" not in df.columns:
        df.reset_index(inplace=True)
    df.set_index(
        ["read_type_id", "chromosome_pk", "barcode_id", "mapping_start_bin_start"],
        inplace=True,
    )
    # df.sort_index(inplace=True)
    df["bin_change"] = df.groupby(
        ["read_type_id", "chromosome_pk", "barcode_id", "mapping_start_bin_start"]
    ).size()
    df_mapping_start_bins = df.loc[~df.index.duplicated()][
        ["bin_change", "run_id", "is_pass"]
    ]
    df_mapping_start_bins["is_start"] = True
    df_mapping_start_bins.index = df_mapping_start_bins.index.rename(
        ["read_type_id", "chromosome_pk", "barcode_id", "bin_position_start"]
    )
    df.reset_index(inplace=True)
    df.set_index(
        ["read_type_id", "chromosome_pk", "barcode_id", "mapping_end_bin_start"],
        inplace=True,
    )
    # df.sort_index(inplace=True)
    df["end_bin_change"] = df.groupby(
        ["read_type_id", "chromosome_pk", "barcode_id", "mapping_end_bin_start"]
    ).size()
    # end_bins are
    df_mapping_end_bins = df.loc[~df.index.duplicated()][
        ["end_bin_change", "run_id", "is_pass"]
    ]
    df_mapping_end_bins["is_start"] = False
    df_mapping_end_bins = df_mapping_end_bins.rename(
        columns={"end_bin_change": "bin_change"}
    )
    df_mapping_end_bins["bin_change"] = -df_mapping_end_bins["bin_change"]
    # If a read ends in a bin, it still covers it, so shift the neg incr to the next bin
    # df_mapping_end_bins["bin_change"] = df_mapping_end_bins["bin_change"].shift().fillna(0)
    df_mapping_end_bins.index = df_mapping_end_bins.index.rename(
        ["read_type_id", "chromosome_pk", "barcode_id", "bin_position_start"]
    )
    results = df_mapping_start_bins.append(df_mapping_end_bins).sort_index()
    # Collapse duplicates on their start + end value
    results["bin_change"] = results.groupby(results.index.names)["bin_change"].transform(np.sum)
    results = results.loc[~results.index.duplicated()]

    results["zero_change"] = results["bin_change"].values.cumsum()
    results.reset_index(inplace=True)
    results["bin_shift"] = results["bin_position_start"].shift(-1)
    results["next_is_start"] = results["is_start"].shift(-1)
    # results = results[results["zero_change"] != 0]

    go_fast2(np.array([[0, 10, 0, 0]], dtype=np.uint32), np.uint32(10))
    results.reset_index(inplace=True)
    positions = results.set_index(["read_type_id", "chromosome_pk"]).groupby(["read_type_id", "chromosome_pk"]).apply(
        lambda x: go_fast2(x[["bin_position_start", "bin_shift", "is_start", "next_is_start"]].fillna(0).values.astype(np.uint32), np.uint32(bin_width)))
    results.set_index(
        ["read_type_id", "chromosome_pk", "barcode_id", "bin_position_start"], inplace=True
    )
    logger.debug(f"Positions {positions.head()}")
    logger.debug(f"positions length {positions.shape}")
    logger.debug(f"positions thing {positions[0]}")
    results = results.loc[~results.index.duplicated()]
    df_old_coverage = fetch_rough_cov_to_update(job_master_id, positions)
    df_to_process = results.append(df_old_coverage, sort=True).sort_index()
    logger.debug("our old coverage")
    logger.debug(df_old_coverage.shape)
    logger.debug(df_old_coverage.head())
    del df_old_coverage
    if "bin_coverage" not in df_to_process.keys():
        df_to_process["bin_coverage"] = 0
    else:
        df_to_process["bin_coverage"] = df_to_process["bin_coverage"].fillna(method="ffill").fillna(0)
    df_to_process["bin_change"] = df_to_process["bin_change"].fillna(0)

    # df_to_process = df_to_process[~df_to_process["reference_id"].isnull()]
    # TODO we don't need this? Can do before not on group by
    df_to_process["summed_bin_change"] = df_to_process.groupby(df_mapping_start_bins.index.names[:-1])[
        "bin_change"
    ].agg(np.cumsum)
    # Don't save where's
    df_to_process = df_to_process[df_to_process["summed_bin_change"] != 0]
    df_to_process["bin_coverage"] = df_to_process["summed_bin_change"] + df_to_process["bin_coverage"]
    df_to_process = df_to_process.loc[~df_to_process.index.duplicated(keep="last")]
    if "ORM" in df_to_process.keys():
        logger.debug("Coverage needs updating")
        df_to_update = df_to_process[~df_to_process["ORM"].isnull()]
        df_to_create = df_to_process[df_to_process["ORM"].isnull()]
        df_to_update["ORM"] = df_to_update.apply(update_orm_objects, axis=1)
        PafRoughCov.objects.bulk_update(
            df_to_update["ORM"].values, ["bin_coverage"], batch_size=5000
        )
    else:
        df_to_create = df_to_process
    if not df_to_create.empty:
        df_to_create.reset_index(inplace=True)
        df_to_create["model_objects"] = df_to_create.apply(
            generate_paf_rough_cov_objects,
            axis=1,
            args=(flowcell_id, job_master_id, reference_id),
        )
        PafRoughCov.objects.bulk_create(
            df_to_create["model_objects"].values.tolist(), batch_size=5000
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
    # The JobMaster object
    job_master = JobMaster.objects.get(pk=job_master_id)
    # Unpack metadata about this job
    last_read_pk = int(job_master.last_read)
    reference_pk = int(job_master.reference.id)
    # We're starting fresh boys
    if last_read_pk is None:
        last_read_pk = 0
    # The debugrmation about the reference for this alignment
    reference_info = ReferenceInfo.objects.get(pk=reference_pk)
    # A dictionary with all the results for each chromosome
    chromosome_dict = {}
    # Get all the chromosomes for this reference
    chromosomes = reference_info.referencelines.all()
    longest_chromosome = chromosomes.aggregate(Max("chromosome_length"))["chromosome_length__max"]
    # For each Chromosome in the chromosomes, set the name of the chromosome as key, and the object as the value
    for chromosome in chromosomes:
        chromosome_dict[chromosome.line_name] = chromosome.id
    # # This gives us a json string that we can query against the alignment server.
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
    if paf_df.empty:
        logger.debug("No mapping results.")
        return
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
        paf_rough_coverage_calculations(df, job_master, longest_chromosome)
    else:
        paf_summary_calculations(df, job_master_id)
        df.reset_index(inplace=True)
        df2 = df
        paf_rough_coverage_calculations(df2, job_master, longest_chromosome)
        if "barcode_name" not in df.keys():
            df.reset_index(inplace=True)
        df["barcode_name"] = "All reads"
        # save results for all reads
        paf_summary_calculations(df, job_master_id)
        paf_rough_coverage_calculations(df2, job_master, longest_chromosome)

    return last_read_pk
