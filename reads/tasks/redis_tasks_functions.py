"""
Contains code relating to tasks that interact with our redis cache for storing reads or aggregations about reads and
 flowcells.
"""
import json
from datetime import datetime

import numpy as np
import pandas as pd
import pytz
import redis
from celery import task
from celery.utils.log import get_task_logger

from artic.task_artic_alignment import run_artic_pipeline
from reads.models import (
    FastqRead,
    Flowcell,
    FlowcellChannelSummary,
    FlowcellStatisticBarcode,
    FlowcellHistogramSummary,
    FlowcellSummaryBarcode,
    JobMaster,
    Run,
    JobType,
    FastqReadType,
    Barcode,
    HistogramSummary,
)

logger = get_task_logger(__name__)
redis_instance = redis.StrictRedis(
    host="127.0.0.1", port=6379, db=0, decode_responses=True
)


def split_flowcell(
    existing_or_new_flowcell, from_flowcell_id, to_flowcell_id, to_flowcell_name, run_id
):
    """
    Split a run from a flowcell and move it and it's reads to a new flowcell
    Parameters
    ----------
    existing_or_new_flowcell: str
        one of new or existing, if the run is being moved to a new or existing flowcell
    from_flowcell_id: int
        PK of the flowcell that the run is being moved from.
    to_flowcell_id: int
        PK of the flowcell that the run is being moved, to, if it exists
    to_flowcell_name: str
        Name of flowcell to create to move run to, if not existing
    run_id: int
        The PK of the run that is being moved between flowcells.

    Returns
    -------

    """
    run = Run.objects.get(pk=run_id)

    if existing_or_new_flowcell == "new":
        from_flowcell = Flowcell.objects.get(pk=from_flowcell_id)
        to_flowcell = Flowcell.objects.create(
            name=to_flowcell_name,
            sample_name=to_flowcell_name,
            owner=from_flowcell.owner,
        )

        JobMaster.objects.create(
            flowcell=to_flowcell,
            job_type=JobType.objects.filter(name="UpdateFlowcellDetails")[0],
            last_read=0,
        )

        logger.info(
            "Moving run {} from flowcell {} to flowcell {}".format(
                run.id, from_flowcell.id, to_flowcell.id
            )
        )

        run.flowcell = to_flowcell
        run.save()

    else:
        from_flowcell = Flowcell.objects.get(pk=from_flowcell_id)
        to_flowcell = Flowcell.objects.get(pk=to_flowcell_id)

        logger.info(
            "Moving run {} from flowcell {} to flowcell {}".format(
                run.id, from_flowcell.id, to_flowcell.id
            )
        )

        run.flowcell = to_flowcell
        run.save()

    move_reads_to_flowcell.delay(run.id, to_flowcell.id, from_flowcell.id)

    return run, from_flowcell, to_flowcell


@task()
def move_reads_to_flowcell(run_id, flowcell_id, from_flowcell_id):
    """
    Move reads from a flowcell run to a different flowcell
    Parameters
    ----------
    run_id: int
        The primary key of the run that generated the reads.
    flowcell_id: int
        The primary key of the flowcell record the reads are being moved to.
    from_flowcell_id: int
        The primary key of the flowcell record the reads are being moved from .

    Returns
    -------
    None
    """
    flowcell = Flowcell.objects.get(pk=flowcell_id)
    old_flowcell = Flowcell.objects.get(pk=from_flowcell_id)
    run = Run.objects.get(pk=run_id)
    chunk_size = 100000
    reads = FastqRead.objects.filter(flowell=old_flowcell, run=run)[0:chunk_size]
    reads_len = reads.count()

    logger.info(
        f"Moving {reads_len} reads from run {run.runid} to flowcell {flowcell.name}"
    )

    if reads:
        # bulk update
        reads.update(flowcell=flowcell)
        # get list of dict of the flowcell
        reads = reads.values()
        update_flowcell.delay(reads)
        logger.info(
            f"Finished moving reads from run {run.runid} to flowcell {flowcell.name}"
        )

    # Update tasks on old flowcell
    if not reads:
        reset_flowcell_info.delay(old_flowcell.id)

    # If there are reads left, call the task again.
    if FastqRead.objects.filter(flowell=old_flowcell, run=run).last():
        move_reads_to_flowcell.apply_async(args=(run_id, flowcell_id, from_flowcell_id))


@task()
def reset_flowcell_info(old_flowcell_id):
    """
    Update the old flowcell base-called summaries and flowcell metadata.
    Parameters
    ----------
    old_flowcell_id: int
        The primary key of the flowcell that had reads moved away from it.
    Returns
    -------
    None

    """
    old_flowcell = Flowcell.objects.get(pk=old_flowcell_id)
    # We need to delete the old objects for this flowcell
    FlowcellSummaryBarcode.objects.filter(flowcell=old_flowcell).delete()
    FlowcellStatisticBarcode.objects.filter(flowcell=old_flowcell).delete()
    FlowcellHistogramSummary.objects.filter(flowcell=old_flowcell).delete()
    FlowcellChannelSummary.objects.filter(flowcell=old_flowcell).delete()

    old_flowcell.start_time = datetime.datetime.now()
    old_flowcell.max_channel = 0
    old_flowcell.number_reads = 0
    old_flowcell.number_runs = 0
    old_flowcell.number_barcodes = 0
    old_flowcell.average_length = 0
    old_flowcell.total_read_length = 0
    old_flowcell.has_fastq = False
    old_flowcell.save()

    # recalculate old flowcell

    chunk_size = 100000

    total_reads = FastqRead.objects.filter(flowcell=old_flowcell).count()

    steps = int((total_reads - total_reads % chunk_size) / chunk_size + 1)

    for i in range(0, steps):

        modstart = i * chunk_size

        modend = (i + 1) * chunk_size

        reads = FastqRead.objects.filter(flowcell=old_flowcell)[modstart:modend]

        reads_len = reads.count()

        update_flowcell.delay(reads.values())


def save_flowcell_statistic_barcode(row):
    """
    Save flowcell statistic barcode rows into redis, applied to a pandas group by object
    Parameters
    ----------
    row: pd.core.series.Series
        The row of the group by object.
    Returns
    -------
    None
    """

    do_whatever_this_does(row, "flowcellStatisticBarcode")


def save_flowcell_summary_barcode(row):
    """
    Save flowcell summary barcode rows into redis, applied to a pandas group by object
    Parameters
    ----------
    row: pd.core.series.Series
        The row of the group by object.
    Returns
    -------
    None

    """
    do_whatever_this_does(row, "flowcellSummaryBarcode")


# TODO dump to redis, message key and retreive
@task()
def update_flowcell(reads_list):
    """
    Update the flowcell metadata, including read counts, yield,max_channel, base-called data summaries
    Parameters
    ----------
    reads_list: list of dict
        List of read dictionaries

    Returns
    -------
    None
    """
    read_df = pd.DataFrame(reads_list)
    read_df["channel"] = read_df["channel"].astype(int)
    read_df["sequence_length"] = read_df["sequence_length"].astype(int)
    read_df = read_df.rename(
        columns={
            "rejected_barcode": "rejected_barcode_id",
            "barcode": "barcode_id",
            "type": "type_id",
            "fastqfile": "fastqfile_id",
        }
    )

    read_df = read_df.drop(columns=["sequence", "quality"])

    results = read_df.groupby(["flowcell_id"]).agg(
        {"read_id": "count", "channel": "max", "sequence_length": "sum"}
    )
    results.reset_index().apply(lambda row: update_flowcell_counts(row), axis=1)

    # We can assume that we will always have only a single flowcell in a b
    # batch of reads to process due to the way that minFQ behaves.
    # Todo: catch problems with multiple flowcells in a single transaction

    read_df["status"] = np.where(read_df["is_pass"] == False, "Fail", "Pass")

    read_df["start_time_round"] = np.array(
        ### Converting to store data in 10 minute windows.
        read_df["start_time"],
        dtype="datetime64[m]",
    )

    read_df["start_time_truncate"] = read_df["start_time_round"].apply(
        lambda dt: datetime(dt.year, dt.month, dt.day, dt.hour, 10 * (dt.minute // 10))
    )

    read_df_all_reads = read_df.copy()
    read_df_all_reads["barcode_name"] = "All reads"

    fastq_df = read_df.append(read_df_all_reads, sort=False)
    fastq_df["barcode_name"] = fastq_df["barcode_name"].fillna("No barcode")
    #
    # Calculates statistics for flowcellSummaryBarcode
    #
    fastq_df_result = fastq_df.groupby(
        ["flowcell_id", "barcode_name", "rejected_barcode_id", "type_id", "is_pass",]
    ).agg(
        {
            "sequence_length": ["min", "max", "sum", "count"],
            "quality_average": ["sum"],
            "channel": ["unique"],
        }
    )

    fastq_df_result.reset_index().apply(
        lambda row: save_flowcell_summary_barcode(row), axis=1
    )

    #
    # Calculates statistics for FlowcellStatisticsBarcode
    #

    fastq_df_result = fastq_df.groupby(
        [
            "flowcell_id",
            "start_time_truncate",
            "barcode_name",
            "type_id",
            "is_pass",
            "rejected_barcode_id",
        ]
    ).agg(
        {
            "sequence_length": ["min", "max", "sum", "count"],
            "quality_average": ["sum"],
            "channel": ["unique"],
        }
    )

    fastq_df_result.reset_index().apply(
        lambda row: save_flowcell_statistic_barcode(row), axis=1,
    )

    fastq_df["bin_index"] = (
        fastq_df["sequence_length"]
        - fastq_df["sequence_length"] % HistogramSummary.BIN_WIDTH
    ) / HistogramSummary.BIN_WIDTH

    fastq_df_result = fastq_df.groupby(
        [
            "flowcell_id",
            "barcode_name",
            "type_id",
            "is_pass",
            "bin_index",
            "rejected_barcode_id",
        ]
    ).agg({"sequence_length": ["sum", "count"]})

    fastq_df_result.reset_index().apply(
        lambda row: save_flowcell_histogram_summary(row), axis=1,
    )
    #
    # Calculates statistics for ChannelSummary
    #
    fastq_df_result = fastq_df.groupby(["flowcell_id", "channel"]).agg(
        {"sequence_length": ["sum", "count"]}
    )
    fastq_df_result.reset_index().apply(
        lambda row: save_flowcell_channel_summary(row), axis=1
    )


def check_if_flowcell_has_streamable_tasks(flowcell_pk):
    """
    Check if the flowcell has streamable tasks.
    Parameters
    ----------
    flowcell_pk: int
        Primary key of the flowcell to check tasks from.
    Returns
    -------
    tasks: list of int
        List of streamable job_masters IDs
    """

    streamable_tasks_pks = [16]
    tasks = JobMaster.objects.filter(
        flowcell_id=flowcell_pk, job_type__id__in=streamable_tasks_pks
    ).values("id", "job_type_id")
    return tasks


def sort_reads_by_flowcell_fire_tasks(reads):
    """
    Sort group reads by flowcell, check if
    Parameters
    ----------
    reads: list of dict
        Reads dictionaries.
    Returns
    -------

    """
    read_df = pd.DataFrame(reads)
    flowcell_gb = read_df.groupby("flowcell_id")
    flowcell_ids = np.unique(read_df["flowcell_id"].values).tolist()

    # TODO talk this through with Matt, what if we only have like 100 reads in a flowcell?
    for flowcell_id in flowcell_ids:
        task_lookups = check_if_flowcell_has_streamable_tasks(flowcell_id)
        flowcell_reads = flowcell_gb.get_group(flowcell_id).to_dict(orient="records")
        for task in task_lookups:
            if task["job_type_id"] == 16:

                run_artic_pipeline.delay(task["id"], flowcell_reads)



@task()
def harvest_reads():
    """
    Celery task to pull all reads JSON accrued out of redis cache and run base-called metadata calculations on it.
    Parameters
    ----------

    Returns
    -------
    None
    """
    # Check that we aren't currently running this task
    if (
        not redis_instance.get("harvesting")
        or not redis_instance.get("harvesting") == "1"
    ):
        redis_instance.set("harvesting", 1)
        reads = list()
        count = redis_instance.scard("reads")

        while count > 0:
            read_chunk = redis_instance.spop("reads")
            read_chunk = json.loads(read_chunk)
            reads.extend(read_chunk)
            count -= 1

        if reads:
            logger.info("Analysing Chunk")
            logger.info("Harvest Reads Updating Flowcell.")
            update_flowcell.delay(reads)
            sort_reads_by_flowcell_fire_tasks(reads)
        redis_instance.set("harvesting", 0)


@task()
def save_reads_bulk(reads):
    """
    Save reads into redis after they arrive from minFQ, and to the database for tasks
     to be run on them later.
    Parameters
    ----------
    reads: list of dict
        A list of reads in dictionary form sent from minFQ
    Returns
    -------
    None

    """

    reads_list = []
    run_dict = {}
    for read in reads:
        run_pk = read.get("run", -1)
        if run_pk not in run_dict and run_pk != -1:
            run = Run.objects.get(pk=run_pk)
            run_dict[run_pk] = run
            read["run_id"] = run.id
            read["flowcell_id"] = run.flowcell.id
        else:
            read["run_id"] = run_dict[run_pk].id
            read["flowcell_id"] = run_dict[run_pk].flowcell.id
        fastq_read = FastqRead(
            read_id=read["read_id"],
            read=read["read"],
            channel=read["channel"],
            barcode_id=read["barcode"],
            rejected_barcode_id=read["rejected_barcode"],
            barcode_name=read["barcode_name"],
            sequence_length=read["sequence_length"],
            quality_average=read["quality_average"],
            sequence=read["sequence"],
            quality=read["quality"],
            is_pass=read["is_pass"],
            start_time=read["start_time"],
            run_id=read["run_id"],
            flowcell_id=read["flowcell_id"],
            type_id=read["type"],
            fastqfile_id=read["fastqfile"],
        )
        reads_list.append(fastq_read)

    # Save reads to redis for later processing of base-called data summaries.
    reads_as_json = json.dumps(reads)
    redis_instance.sadd("reads", reads_as_json)
    # Bulk create the entries
    FastqRead.objects.bulk_create(reads_list, batch_size=1000)


def get_values_and_delete_key_redis_hash(r, key):
    """
    Safely delete a redis hash (Dictionary)
    Parameters
    ----------
    r: redis instance
    key: str
        key to fetch and delete safely.

    Returns
    -------
    key value
    """
    p = r.pipeline()
    p.multi()
    p.hgetall(key)
    p.delete(key)
    returnvalue = p.execute()[0]
    return returnvalue


def update_flowcell_counts(row):
    """
    Update flowcell data counts, has_fastq, number_reads, yield, max channel
    Parameters
    ----------
    row: pandas.core.series.Series
        A group by object row.
    Returns
    -------
    None

    """
    flowcell_id = int(row["flowcell_id"])
    # TODO this is where we would set no fastq of we WEREN't saving te data but were uploading it
    redis_instance.set(f"{flowcell_id}_has_fastq", 1)
    # Increment by the count of the read_ids on this flowcell for total read number
    redis_instance.incrby(f"{flowcell_id}_number_reads", int(row["read_id"]))
    #  Increment by the sum of the sequence length to calculate yield
    redis_instance.incrby(
        f"{flowcell_id}_total_read_length", int(row["sequence_length"]),
    )
    # Current largest channel
    saved_max_channel = redis_instance.get(f"{flowcell_id}_max_channel")
    if saved_max_channel:
        # Update to new, larger channel
        if (row["channel"]) > int(saved_max_channel):
            redis_instance.set(f"{flowcell_id}_max_channel", int(row["channel"]))
    else:
        # Set initial max channel
        redis_instance.set(f"{flowcell_id}_max_channel", int(row["channel"]))


def save_flowcell_histogram_summary(row):
    """
    Save aggregated flowcell histogram summary data into the redis cache.
    Parameters
    ----------
    row: pandas.core.series.Series
        The row of a pandas group by object

    Returns
    -------
    None

    """
    flowcell_id = int(row["flowcell_id"][0])
    # flowcell = Flowcell.objects.get(pk=flowcell_id)
    barcode_name = row["barcode_name"][0]
    read_type_name = FastqReadType.objects.get(pk=row["type_id"][0]).name
    # read_type_name = row["type__name"][0]
    rejection_status = Barcode.objects.get(pk=row["rejected_barcode_id"][0]).name
    # rejection_status = row["rejected_barcode__name"][0]
    status = row["is_pass"][0]
    bin_index = int(row["bin_index"][0])
    sequence_length_sum = int(row["sequence_length"]["sum"])
    read_count = int(row["sequence_length"]["count"])

    unique_key = "{}_flowcellHistogramSummary_{}_{}_{}_{}_{}".format(
        flowcell_id, barcode_name, rejection_status, read_type_name, status, bin_index
    )
    p = redis_instance.pipeline()
    p.multi()
    p.hset(unique_key, "bin_index", bin_index)
    p.hset(unique_key, "barcode_name", barcode_name)
    p.hset(unique_key, "rejection_status", rejection_status)
    p.hset(unique_key, "read_type_name", read_type_name)
    p.hset(unique_key, "status", str(status))
    p.hincrby(unique_key, "read_length", sequence_length_sum)
    p.hincrby(unique_key, "read_count", read_count)
    p.execute()


def save_flowcell_channel_summary(row):
    """
    Increment a channel presence count for each channel in the redis cache.
    Parameters
    ----------
    row: pandas.core.series.Series

    Returns
    -------
    None
    """
    flowcell_id = row["flowcell_id"][0]
    channel_number = row["channel"][0]
    channel_read_count = int(row["sequence_length"]["count"])
    channel_yield = int(row["sequence_length"]["sum"])
    redis_instance.incrby(
        f"{flowcell_id}_{channel_number}_read_count", channel_read_count,
    )
    redis_instance.incrby(
        f"{flowcell_id}_{channel_number}_pore_yield", channel_yield,
    )


def do_whatever_this_does(row, unique_key_base):
    """
    Store the aggregated results of the flowcell statistic calculations and flowcell summary barcode
     calculations in redis.
    Parameters
    ----------
    row: pandas.core.series.Series
        A group by series
    unique_key_base: str
        Base for the unique redis key string. One of flowcellStatisticBarcode or flowcellSummaryBarcode

    Returns
    -------

    """
    flowcell_id = int(row["flowcell_id"][0])
    type_name = FastqReadType.objects.get(pk=row["type_id"][0]).name
    barcode_name = row["barcode_name"][0]
    rejection_status = Barcode.objects.get(pk=row["rejected_barcode_id"][0]).name
    status = row["is_pass"][0]
    sequence_length_sum = int(row["sequence_length"]["sum"])
    sequence_length_max = int(row["sequence_length"]["max"])
    sequence_length_min = int(row["sequence_length"]["min"])
    quality_average_sum = int(row["quality_average"]["sum"])
    read_count = int(row["sequence_length"]["count"])
    channels = row["channel"]["unique"]
    if unique_key_base == "flowcellStatisticBarcode":
        utc = pytz.utc
        start_time = utc.localize(row["start_time_truncate"][0])
        unique_key = (
            f"{flowcell_id}_{unique_key_base}_{start_time}_{barcode_name}"
            f"_{rejection_status}_{type_name}_{status}"
        )
    else:
        unique_key = f"{flowcell_id}_{unique_key_base}_{barcode_name}_{rejection_status}_{type_name}_{status}"
    prev_channel_list = redis_instance.hget(unique_key, "channel_presence")
    min_prev = redis_instance.hget(unique_key, "min_length")
    max_prev = redis_instance.hget(unique_key, "max_length")
    p = redis_instance.pipeline()
    p.multi()
    if unique_key_base == "flowcellStatisticBarcode":
        p.hset(unique_key, "sample_time", str(start_time))
    p.hset(unique_key, "barcode_name", barcode_name)
    p.hset(unique_key, "rejection_status", rejection_status)
    p.hset(unique_key, "type_name", type_name)
    p.hset(unique_key, "status", str(status))
    p.hincrby(unique_key, "total_length", sequence_length_sum)
    p.hincrby(unique_key, "quality_sum", quality_average_sum)
    p.hincrby(unique_key, "read_count", read_count)

    if min_prev:
        min_prev = int(min_prev)
        if min_prev > sequence_length_min:
            p.hset(unique_key, "min_length", sequence_length_min)
    else:
        p.hset(unique_key, "min_length", sequence_length_min)

    if max_prev:
        max_prev = int(max_prev)
        if max_prev < sequence_length_max:
            p.hset(unique_key, "max_length", sequence_length_max)
    else:
        p.hset(unique_key, "max_length", sequence_length_max)

    if prev_channel_list:
        channel_list = list(prev_channel_list)
    else:
        channel_list = list("0" * 3000)

    for c in channels:
        channel_list[int(c) - 1] = "1"

    channel_presence = "".join(channel_list)

    p.hset(unique_key, "channel_presence", channel_presence)
    p.execute()


def get_values_and_delete_redis_key(r, key):
    """
    Safely delete a key
    Parameters
    ----------
    r: redis instance
    key: key to fetch and delete safely.

    Returns
    -------
    key value
    """
    p = r.pipeline()
    p.multi()
    p.get(key)
    p.delete(key)
    returnvalue = p.execute()[0]
    return returnvalue


def scan_keys(r, pattern):
    """
    Returns a list of all the keys matching a given pattern in a redis index.
    Parameters
    ----------
    r: redis.client.Redis
        Instance allowing us to access the redis cache
    pattern: str
        Pattern string to match against keys for.

    Returns
    -------

    """
    result = []
    cur, keys = r.scan(cursor=0, match=pattern, count=2)
    result.extend(keys)
    while cur != 0:
        cur, keys = r.scan(cursor=cur, match=pattern, count=2)
        result.extend(keys)
    return result
