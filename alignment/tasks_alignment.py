from __future__ import absolute_import, unicode_literals

from collections import namedtuple

import numpy as np
import pandas as pd
import redis
from celery.utils.log import get_task_logger
from django.db.models import F, Max
from numba import njit

from alignment.mapper import MAP
from alignment.models import PafRoughCov, PafSummaryCov, PafRoughCovIntermediate
from minotourapp.celery import app
from reads.models import FastqRead, JobMaster, Barcode
from reference.models import ReferenceInfo

# Set up the logger to write to logging file
logger = get_task_logger(__name__)
redis_instance = redis.StrictRedis(
    host="127.0.0.1", port=6379, db=0, decode_responses=True
)


@app.task
def aggregate_intermediate_table():
    """
    Aggregate the intermediate rough cov bins and save them
    Returns
    -------

    """
    fetch_intmd_rough_cov_to_aggregate()


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
        print(f"Fetching reads in chunks of {chunk_size} for alignment.")
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
            print("No fastq found. Skipping this iteration.")
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

    print(
        f"Flowcell id: {flowcell.id} - Running minimap2 on flowcell {flowcell.name}"
    )
    print(f"Flowcell id: {flowcell.id} - job_master_id {job_master.id}")
    print(f"Flowcell id: {flowcell.id} - reference {job_master.reference.name}")
    print(f"Flowcell id: {flowcell.id} - last_read {job_master.last_read}")
    print(f"Flowcell id: {flowcell.id} - read_count {job_master.read_count}")
    print(f"Flowcell id: {flowcell.id} - number of reads found {read_count}")
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


def save_summary_coverage(row, job_master_id, reference_id, reference_name):
    """
    Save the paf store summaries
    Parameters
    ----------
    row: pd.core.series.Series
        The row of the dataframe we are saving
    job_master_id: int
        The primary key of the databases task record
    reference_id: int
        The primary key of the reference record
    reference_name: str
        The name of the reference
    Returns
    -------

    """
    paf_summary_cov, created = PafSummaryCov.objects.get_or_create(
        job_master_id=job_master_id,
        barcode_name=row["barcode_name"],
        barcode_id=row["barcode_id"],
        chromosome_id=row["chromosome_pk"],
        reference_line_length=row["chromosome_length"],
        read_type_id=row["read_type_id"],
        reference_pk=reference_id,
        chromosome_pk=row["chromosome_pk"],
        chromosome_name=row["tsn"],
        reference_name=reference_name,
    )
    paf_summary_cov.total_yield += row["yield"]
    paf_summary_cov.read_count += row["read_count"]
    paf_summary_cov.coverage = round(
        paf_summary_cov.total_yield / int(row["chromosome_length"]), 2
    )
    paf_summary_cov.average_read_length = round(
        paf_summary_cov.total_yield / paf_summary_cov.read_count, 2
    )
    paf_summary_cov.save()


def update_final_rough_coverage(df_to_process):
    """
    Update the final rough coverage of by combining existing values with new values and saving them

    Parameters
    ----------
    df_to_process: pandas.core.frame.DataFrame
        The datframe of Intermediate PafRoughCovValues
    Returns
    -------

    """
    print("update final rough cov")
    if "bin_coverage" not in df_to_process.keys():
        df_to_process["bin_coverage"] = 0
    else:
        df_to_process["bin_coverage"] = (
            df_to_process["bin_coverage"].fillna(method="ffill").fillna(0)
        )
    df_to_process["bin_change"] = df_to_process["bin_change"].fillna(0)
    print(df_to_process.head())
    print(df_to_process.index)
    print(df_to_process.keys())
    df_to_process["summed_bin_change"] = df_to_process.groupby(
        df_to_process.index.names[:-2]
    )["bin_change"].agg(np.cumsum)
    df_to_process["bin_coverage"] = (
        df_to_process["summed_bin_change"] + df_to_process["bin_coverage"]
    )
    df_to_process = df_to_process.loc[~df_to_process.index.duplicated(keep="last")]
    print(df_to_process.head())
    if "ORM" in df_to_process.keys():
        print("Coverage needs updating")
        df_to_update = df_to_process[~df_to_process["ORM"].isnull()]
        df_to_create = df_to_process[df_to_process["ORM"].isnull()]
        # Don't save where the change is 0
        df_to_update = df_to_update[df_to_update["summed_bin_change"].ne(0)]
        df_to_update["ORM"] = df_to_update.apply(update_orm_objects, axis=1)
        PafRoughCov.objects.bulk_update(
            df_to_update["ORM"].values, ["bin_coverage"], batch_size=5000
        )
    else:
        df_to_create = df_to_process
    print(df_to_create)
    print(df_to_create.keys())
    if not df_to_create.empty:
        print("creating")
        df_to_create.reset_index(inplace=True)
        df_to_create["model_objects"] = df_to_create.apply(
            generate_paf_rough_cov_objects,
            axis=1,
            args=(PafRoughCov, "chromosome_name"),
        )
        PafRoughCov.objects.bulk_create(
            df_to_create["model_objects"].values.tolist(), batch_size=5000
        )

def fetch_final_rough_cov_to_update(df_aggregate_intmd):
    """
    Fetch the final data from the paf rough cov table to update
    Parameters
    ----------
    df_aggregate_intmd: pd.core.frame.DataFrame
        The aggregated dataframe of all the results in the intermediate table, grouped by job_master, position,
        chromosome and read type
    Returns
    -------

    """
    go_fast2(np.array([[0, 10, 0, 0]], dtype=np.uint32), np.uint32(10))
    df_aggregate_intmd.reset_index(inplace=True)
    bin_width = 10
    positions = (
        df_aggregate_intmd.set_index(["job_master_id", "read_type_id", "chromosome_pk"])
        .groupby(["job_master_id", "read_type_id", "chromosome_pk"])
        .apply(
            lambda x: go_fast2(
                x[["bin_position_start", "bin_shift", "is_start", "next_is_start"]]
                .fillna(0)
                .values.astype(np.uint32),
                np.uint32(bin_width),
            )
        )
    )
    print(f"Positions {positions.head()}")
    print(f"positions length {positions.shape}")
    print(f"positions thing {positions[0]}")
    querysets = (
        positions.groupby(["job_master_id", "chromosome_pk"]).apply(lambda x: PafRoughCov.objects.filter(
            chromosome_id=x.name[1], job_master_id=x.name[0], bin_position_start__in=x.values[0]
        )).values
    )
    if not querysets[0]:
        print("No previous bins found")
        df_old_coverage = pd.DataFrame()
        df_aggregate_intmd.set_index([
            "read_type_id",
            "chromosome_pk",
            "barcode_id",
            "bin_position_start",
            "job_master_id"
        ], inplace=True)
        df_to_process = df_aggregate_intmd.append(df_old_coverage, sort=True).sort_index()
        update_final_rough_coverage(df_to_process)
        return

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
            columns=["_state", "grouprun_id", "bin_position_end", "is_pass", "chromosome_id"], inplace=True
        )
        # df_old_coverage.rename(columns={"chromosome_id": "chromosome_pk"}, inplace=True)
        df_old_coverage["ORM"] = df_old_coverage["id"].map(lookup_dict)
        df_old_coverage.set_index(
            [
                "read_type_id",
                "chromosome_pk",
                "barcode_id",
                "bin_position_start",
                "job_master_id"
            ], inplace=True
        )
        print(df_old_coverage)
        df_old_coverage = df_old_coverage.loc[~df_old_coverage.index.duplicated()]
    df_aggregate_intmd.set_index([
                "read_type_id",
                "chromosome_pk",
                "barcode_id",
                "bin_position_start",
                "job_master_id"
            ], inplace=True)
    df_to_process = df_aggregate_intmd.append(df_old_coverage, sort=True).sort_index()
    update_final_rough_coverage(df_to_process)


def fetch_intmd_rough_cov_to_aggregate():
    """
    Returns a DataFrame of PafRoughCov bins we have data for
    Returns
    -------
    pd.core.frame.DataFrame
    """
    querysets = PafRoughCovIntermediate.objects.all()
    if not querysets:
        print("No Intermediate data found, exiting.")
        return
    # lookup_dict = {prc.id: prc for prc in querysets}
    df_to_aggregate = pd.DataFrame.from_records(querysets.values())
    querysets.delete()
    # df_to_aggregate["ORM"] = df_to_aggregate["id"].map(lookup_dict)
    df_to_aggregate["bin_shift"] = df_to_aggregate["bin_position_start"].shift(-1)
    df_to_aggregate["next_is_start"] = df_to_aggregate["is_start"].shift(-1)
    df_to_aggregate = df_to_aggregate.set_index(
        [
            "job_master_id",
            "read_type_id",
            "chromosome_pk",
            "barcode_id",
            "bin_position_start",
        ]
    )
    fetch_final_rough_cov_to_update(df_to_aggregate)

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


def generate_paf_rough_cov_objects(
    row, model, key="tsn"
):
    """
    Generate Paf summary coverage ORM objects on a dataframe row
    Parameters
    ----------
    row: pd.core.series.Series
        The row of a dataframe
    model: PafRoughCovIntermediate or PafroughCov
        The model we are updating or saving
    key: str
        The lookup key for the chromosome name
    Returns
    -------
    alignment.models.PafSummaryCov

    """#
    a = model(
        job_master_id=row["job_master_id"],
        run_id=row["run_id"],
        flowcell_id=row["flowcell_id"],
        read_type_id=row["read_type_id"],
        is_pass=row["is_pass"],
        barcode_id=row["barcode_id"],
        reference_id=row["reference_id"],
        chromosome_id=row["chromosome_pk"],
        bin_position_start=row["bin_position_start"],
        chromosome_pk=row["chromosome_pk"],
        reference_pk=row["reference_id"],
        reference_name=row["reference_name"],
        chromosome_length=row["chromosome_length"],
        chromosome_name=row[key],
    )
    if model == PafRoughCovIntermediate:
        a.bin_change = row["bin_change"]
        a.is_start = row["is_start"]

    else:
        a.bin_coverage = row["bin_coverage"]
    return a



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
    print("rough coverage calculations")
    flowcell_id = int(job_master.flowcell.id)
    reference_id = int(job_master.reference.id)
    reference_name = str(job_master.reference.name)
    job_master_id = int(job_master.id)
    bin_width = 10
    # if longest_chromosome / bin_width > 1000000:
    #     bin_width = round(longest_chromosome / 1000000 / 10) * 10
    # if bin_width < 10:
    #     bin_width = 10
    rounder = bin_width / 2
    df["mapping_start_bin_start"] = (
        df["mapping_start"].sub(rounder).div(bin_width).round().mul(bin_width)
    )
    df["mapping_end_bin_start"] = (
        df["mapping_end"].sub(rounder).div(bin_width).round().mul(bin_width)
    )
    if len(df.index.names) > 1:
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
        ["bin_change", "run_id", "is_pass", "chromosome_length", "tsn"]
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
        ["end_bin_change", "run_id", "is_pass", "chromosome_length", "tsn"]
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
    results["bin_change"] = results.groupby(results.index.names)[
        "bin_change"
    ].transform(np.sum)
    print(results.head())
    results = results.loc[~results.index.duplicated()]
    results.reset_index(inplace=True)
    # compile the code that runs this, so it will run super fast for the next one
    results.set_index(
        ["read_type_id", "chromosome_pk", "barcode_id", "bin_position_start"],
        inplace=True,
    )
    results = results.loc[~results.index.duplicated()]
    results["flowcell_id"] = flowcell_id
    results["job_master_id"] = job_master_id
    results["reference_name"] = reference_name
    results["reference_id"] = reference_id
    print("results")
    if not results.empty:
        results.reset_index(inplace=True)
        results["model_objects"] = results.apply(
            generate_paf_rough_cov_objects,
            axis=1,
            args=(PafRoughCovIntermediate,),
        )
        PafRoughCovIntermediate.objects.bulk_create(
            results["model_objects"].values.tolist(), batch_size=5000
        )


def paf_summary_calculations(df, job_master_id, reference):
    """
    Do dataframe aggregations on a per chromosome, barcode and readtype combination level
    and save the resulting paf_summary_cov
    Parameters
    ----------
    df: pd.core.frame.DataFrame
        Merged dataframe of paf output and fastq dataframe
    job_master_id: int
        The primary key of this task database record
    reference: reference.models.ReferenceInfo
        The reference ORM record
    Returns
    -------

    """
    reference_name = str(reference.name)
    reference_pk = int(reference.id)

    if "type__name" in df.keys():
        df.set_index(["type__name", "tsn", "barcode_name"], inplace=True)
    df[["yield", "read_count"]] = df.groupby(
        level=["type__name", "tsn", "barcode_name"]
    ).agg({"mapping_length": np.sum, "sequence_length": np.size})
    df["chromosome_length"] = df["chromosome_length"].astype(int)
    df = df.loc[~df.index.duplicated()]
    df.reset_index(inplace=True)
    print("summary df being saved is:")
    print(df.head())
    print(df.keys())
    df.apply(
        save_summary_coverage,
        args=(job_master_id, reference_pk, reference_name),
        axis=1,
    )


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
    longest_chromosome = chromosomes.aggregate(Max("chromosome_length"))[
        "chromosome_length__max"
    ]
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
        print("No mapping results.")
        return
    del test_df
    # return test_df

    paf_df = paf_df.rename(
        columns={
            0: "read_id",
            1: "sequence_length",
            5: "tsn",
            6: "chromosome_length",
            7: "mapping_start",
            8: "mapping_end",
            11: "mq",
        }
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
    # All reads barcode Id
    run_id = np.unique(df["run_id"].values)[0]
    all_reads_barcode_id = int(Barcode.objects.get(run_id=run_id, name="ALl reads").id)
    # TODO should these go to a chain, final in chain deletes the paf results in redis
    if "No barcode" in fasta_df_barcode["barcode_name"].unique():
        df["barcode_name"] = "All reads"
        df["barcode_id"] = all_reads_barcode_id
        paf_summary_calculations(df, job_master_id, reference_info)
        paf_rough_coverage_calculations(df, job_master, longest_chromosome)
    else:
        paf_summary_calculations(df, job_master_id, reference_info)
        df.reset_index(inplace=True)
        df2 = df
        paf_rough_coverage_calculations(df2, job_master, longest_chromosome)
        if len(df.index.names) > 1:
            df.reset_index(inplace=True)
        df["barcode_name"] = "All reads"
        df["barcode_id"] = all_reads_barcode_id
        print(df.keys())
        print(df.head())
        # save results for all reads
        paf_summary_calculations(df, job_master_id, reference_info)
        paf_rough_coverage_calculations(df2, job_master, longest_chromosome)

    return last_read_pk
