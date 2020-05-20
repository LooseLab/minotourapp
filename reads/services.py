"""
Services to help run some of the tasks in the web app. Why is this not stored there?
"""
import pytz
import datetime
import pandas as pd
import numpy as np
import decimal
import os
import subprocess

from celery.task import task
from celery.utils.log import get_task_logger
from collections import defaultdict
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.db import IntegrityError, transaction

from django.conf import settings

from textwrap import wrap

from reads.models import (
    FlowcellStatisticBarcode,
    FlowcellSummaryBarcode,
    FlowcellHistogramSummary,
    FlowcellChannelSummary,
    HistogramSummary,
    Flowcell,
    Run,
    MinionMessage,
    FastqRead, FastqReadType, Barcode)

from alignment.models import PafRoughCov, PafStore, PafSummaryCov
from reference.models import ReferenceInfo, ReferenceLine

from communication.models import Message, NotificationConditions

import redis

import json

redis_instance = redis.StrictRedis(host="127.0.0.1",
                                  port=6379, db=0,decode_responses=True)


logger = get_task_logger(__name__)

# @task(rate_limit="2/m")

def processreads(reads):
    reads_list = []
    run_dict = {}
    for read in reads:
        run_pk = read.get("run", -1)
        if run_pk not in run_dict and run_pk != -1:
            run = Run.objects.get(pk=run_pk)
            run_dict[run_pk] = run
            read["run"] = run
            read["flowcell"] = run.flowcell
        else:
            read["run"] = run_dict[run_pk]
            read["flowcell"] = run_dict[run_pk].flowcell
        fastqread = FastqRead(
            read_id=read['read_id'],
            read=read['read'],
            channel=read['channel'],
            barcode_id=read['barcode'],
            rejected_barcode_id=read['rejected_barcode'],
            barcode_name=read['barcode_name'],
            sequence_length=read['sequence_length'],
            quality_average=read['quality_average'],
            sequence=read['sequence'],
            quality=read['quality'],
            is_pass=read['is_pass'],
            start_time=read['start_time'],
            run=read["run"],
            flowcell=read["flowcell"],
            type_id=read['type'],
            fastqfile_id=read['fastqfile']
        )
        reads_list.append(fastqread)
    return reads_list


@task()
def save_reads_bulk(reads):
    """This task will save reads into redis for later processing."""
    ### Save reads to redis for later processing.
    reads_as_json = json.dumps(reads)
    redis_instance.sadd("reads", reads_as_json)
    reads_list = processreads(reads)
    FastqRead.objects.bulk_create(reads_list, batch_size=500)



@task(serializer="pickle")
def check_for_jobs(readDF,flowcell):
    """
    This function checks for jobs that can be run on these data and passes them on to tasks as appropriate.
    Parameters
    ----------
    readDF A dataframe containing all the read data
    flowcell The flowcell id to check

    Returns
    -------

    """
    flowcell = Flowcell.objects.get(pk=flowcell)
    for job in flowcell.flowcelljobs.all():
        if job.job_type.name == "Minimap2":
            align_reads.delay(readDF,flowcell)
    pass


@task()
def harvestreads():
    """
    Celery task to bulk create the reads
    Parameters
    ----------
    reads: list of dict
        List of FastqRead dicts to save.

    Returns
    -------
    None
    """
    if not redis_instance.get("harvesting") or not redis_instance.get("harvesting")=="1":
        redis_instance.set("harvesting",1)
        reads = list()
        count = redis_instance.scard("reads")
        logger.info("Harvest Reads had seen {} batches.".format(count))
        while (count > 0):
            read_chunk = redis_instance.spop("reads")
            read_chunk = json.loads(read_chunk)
            reads.append(read_chunk)
            count-=1
        if len(reads)>0:
            reads_list = []
            #run_dict = {}
            for chunk in reads:
                logger.info("Harvesting Chunk")
                reads_list += processreads(chunk)
            logger.info("Harvest Reads Updating Flowcell.")
            if len(reads_list) > 0:
                update_flowcell.delay(reads_list)
        redis_instance.set("harvesting", 0)

@task(serializer="pickle")
def update_flowcell(reads_list):
    """
    Celery task to update the flowcell record for a specific batch of reads.

    We need to do anything that the current flowcell task does based on a query of reads in the database.
    This include:
    set max channel
    update read count

    It will also handle updating some of the chan calc parameters.
    Parameters
    ----------
    readDF

    Returns
    -------

    """
    # This manipulation returns a flowcell_id, read count and max channel for each flowcell in the dataFrame
    #readDF = pd.DataFrame([o.__dict__ for o in reads_list]).drop(['_state'], axis=1)
    readDF = pd.DataFrame([o.__dict__ for o in reads_list])
    readDF['channel'] = readDF['channel'].astype(int)
    readDF['sequence_length'] = readDF['sequence_length'].astype(int)
    ### Lets try mapping some sequence here.
    flowcells = readDF['flowcell_id'].unique()
    for flowcell in flowcells:
        check_for_jobs.delay(readDF[readDF['flowcell_id'] == flowcell], flowcell)

    # This loses all the sequence data etc - so if we want to do anything with mappings etc we want to do it now.
    readDF = readDF.drop(columns=['sequence','quality'])

    results = readDF.groupby(['flowcell_id']).agg({"read_id": "count", "channel": "max","sequence_length":"sum"})
    results.reset_index().apply(
        lambda row: update_flowcell_counts(row), axis=1
    )


    ## We can assume that we will always have only a single flowcell in a batch of reads to process due to the way that minFQ behaves.
    ## Todo: catch problems with multiple flowcells in a single transaction

    readDF["status"] = np.where(
        readDF["is_pass"] == False, "Fail", "Pass"
    )

    readDF["start_time_round"] = np.array(
        ### Converting to store data in 10 minute windows.
        readDF["start_time"], dtype="datetime64[m]"
    )

    readDF["start_time_truncate"] = readDF['start_time_round'].apply(
        lambda dt: datetime.datetime(dt.year, dt.month, dt.day, dt.hour, 10 * (dt.minute // 10)))

    readDF_allreads = readDF.copy()
    readDF_allreads["barcode_name"] = "All reads"

    fastq_df = readDF.append(readDF_allreads,sort=False)
    fastq_df["barcode_name"] = fastq_df["barcode_name"].fillna("No barcode")
    #
    # Calculates statistics for flowcellSummaryBarcode
    #

    fastq_df_result = fastq_df.groupby(
        [
            "flowcell_id",
            "barcode_name",
            "rejected_barcode_id",
            "type_id",
            "is_pass",
        ]
    ).agg(
        {
            "sequence_length": ["min", "max", "sum", "count"],
            "quality_average": ["sum"],
            "channel": ["unique"],
        }
    )

    fastq_df_result.reset_index().apply(
        lambda row: save_flowcell_summary_barcode_async(row), axis=1
    )

    #
    # Calculates statistics for FlowcellStatisticsBarcode
    #
    # todo: implement this.

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
        lambda row: save_flowcell_statistic_barcode_async(row),
        axis=1,
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
        lambda row: save_flowcell_histogram_summary_async(row),
        axis=1,
    )

    #
    # Calculates statistics for ChannelSummary
    #
    fastq_df_result = fastq_df.groupby(["flowcell_id","channel"]).agg(
        {"sequence_length": ["sum", "count"]}
    )
    #print (fastq_df_result)

    fastq_df_result.reset_index().apply(
        lambda row: save_flowcell_channel_summary_async(row), axis=1
    )




#@task(serializer="pickle")
def update_flowcell_counts(row):
    """Playing with redis"""
    redis_instance.set("{}_has_fastq".format(int(row["flowcell_id"])), 1)
    redis_instance.incrby("{}_number_reads".format(int(row["flowcell_id"])),int(row["read_id"]))
    redis_instance.incrby("{}_total_read_length".format(int(row["flowcell_id"])),int(row["sequence_length"]))
    saved_max_channel = redis_instance.get("{}_max_channel".format(int(row["flowcell_id"])))
    if saved_max_channel:
            if (row["channel"])>int(saved_max_channel):
                redis_instance.set("{}_max_channel".format(int(row["flowcell_id"])),int(row["channel"]))
    else:
        redis_instance.set("{}_max_channel".format(int(row["flowcell_id"])),int(row["channel"]))

#@task(serializer="pickle")
def save_flowcell_statistic_barcode_async(row):
    """
    Shifting this to store in redis
    Parameters
    ----------
    row

    Returns
    -------

    """

    flowcell_id = int(row["flowcell_id"][0])
    #flowcell = Flowcell.objects.get(pk=flowcell_id)
    utc = pytz.utc
    barcode_name = row["barcode_name"][0]
    type_name = FastqReadType.objects.get(pk=row["type_id"][0]).name
    start_time = utc.localize(row["start_time_truncate"][0])
    rejection_status = Barcode.objects.get(pk=row["rejected_barcode_id"][0]).name
    status = row["is_pass"][0]
    sequence_length_sum = int(row["sequence_length"]["sum"])
    sequence_length_max = int(row["sequence_length"]["max"])
    sequence_length_min = int(row["sequence_length"]["min"])
    quality_average_sum = int(row["quality_average"]["sum"])
    read_count = int(row["sequence_length"]["count"])
    channels = row["channel"]["unique"]

    #our unique key for this is going to be:
    uniquekey = "{}_flowcellStatisticBarcode_{}_{}_{}_{}_{}".format(flowcell_id,start_time,barcode_name,rejection_status,type_name,status)
    ## Make sure these are all together.
    min_prev = redis_instance.hget(uniquekey, "min_length")
    max_prev = redis_instance.hget(uniquekey, "max_length")
    prev_channel_list = redis_instance.hget(uniquekey, "channel_presence")

    p = redis_instance.pipeline()
    p.multi()
    p.hset(uniquekey, "sample_time", str(start_time))
    p.hset(uniquekey, "barcode_name", barcode_name)
    p.hset(uniquekey, "rejection_status", rejection_status)
    p.hset(uniquekey, "type_name", type_name)
    p.hset(uniquekey, "status", str(status))
    p.hincrby(uniquekey, "total_length", sequence_length_sum)
    p.hincrby(uniquekey, "quality_sum", quality_average_sum)
    p.hincrby(uniquekey, "read_count", read_count)

    if min_prev:
        min_prev = int(min_prev)
        if min_prev > sequence_length_min:
            p.hset(uniquekey, "min_length", sequence_length_min)
    else:
        p.hset(uniquekey, "min_length", sequence_length_min)
    if max_prev:
        max_prev = int(max_prev)
        if max_prev < sequence_length_max:
            p.hset(uniquekey, "max_length", sequence_length_max)
    else:
        p.hset(uniquekey, "max_length", sequence_length_max)
    if prev_channel_list:
        channel_list = list(prev_channel_list)
    else:
        channel_list = list("0"*3000)

    for c in channels:
        channel_list[int(c) - 1] = "1"
    channel_presence = "".join(channel_list)
    p.hset(uniquekey,"channel_presence",channel_presence)
    p.execute()

    """
    flowcellStatisticBarcode, created = FlowcellStatisticBarcode.objects.get_or_create(
        flowcell=flowcell,
        sample_time=start_time,
        barcode_name=barcode_name,
        rejection_status=rejection_status,
        read_type_name=type_name,
        status=status,
    )

    flowcellStatisticBarcode.total_length += sequence_length_sum
    flowcellStatisticBarcode.quality_sum += decimal.Decimal(quality_average_sum)
    flowcellStatisticBarcode.read_count += read_count

    if flowcellStatisticBarcode.max_length < sequence_length_max:
        flowcellStatisticBarcode.max_length = sequence_length_max

    if flowcellStatisticBarcode.min_length < sequence_length_min:
        flowcellStatisticBarcode.min_length = sequence_length_min

    channel_list = list(flowcellStatisticBarcode.channel_presence)

    for c in channels:

        channel_list[int(c) - 1] = "1"

    flowcellStatisticBarcode.channel_presence = "".join(channel_list)
    flowcellStatisticBarcode.channel_count = len(channels)

    flowcellStatisticBarcode.save()
    
    """

#@task(serializer="pickle")
def save_flowcell_channel_summary_async(row):
    """
        Convert to redis.
        Model is store under
        flowcell_id_channel_read_length
        flowcell_id_channel_read_count

        Parameters
        ----------
        row

        Returns
        -------

        """

    redis_instance.incrby("{}_{}_read_count".format(int(row["flowcell_id"]),int(row["channel"][0])), int(row["sequence_length"]["count"]))
    redis_instance.incrby("{}_{}_read_length".format(int(row["flowcell_id"]),int(row["channel"][0])), int(row["sequence_length"]["sum"]))

    """
    Save flowcell channel summary into the database row-wise on a pandas dataframe. Used for the channel visualisations.
    :param flowcell_id: Primary key of the flowcell database entry
    :type flowcell_id: int
    :param row: a pandas dataframe row
    :type row: pandas.core.series.Series
    :return: None
    """
    """
    flowcell_id = int(row["flowcell_id"][0])
    flowcell = Flowcell.objects.get(pk=flowcell_id)
    channel = int(row["channel"][0])
    sequence_length_sum = int(row["sequence_length"]["sum"])
    read_count = int(row["sequence_length"]["count"])
    flowcellChannelSummary, created = FlowcellChannelSummary.objects.get_or_create(
        flowcell=flowcell, channel=channel
    )
    flowcellChannelSummary.read_length += sequence_length_sum
    flowcellChannelSummary.read_count += read_count
    flowcellChannelSummary.save()
    """

#@task(serializer="pickle")
def save_flowcell_histogram_summary_async(row):
    """
    Save flowcell histogram barcode objects into the database row-wise on a pandas dataframe.
    Data for the histograms on the chancalc page.
    :param flowcell_id: Primary key of the flowcell database entry
    :type flowcell_id: int
    :param row: a pandas dataframe row
    :type row: pandas.core.series.Series
    :return: None
    """
    flowcell_id = int(row["flowcell_id"][0])
    #flowcell = Flowcell.objects.get(pk=flowcell_id)
    barcode_name = row["barcode_name"][0]
    read_type_name = FastqReadType.objects.get(pk=row["type_id"][0]).name
    #read_type_name = row["type__name"][0]
    rejection_status = Barcode.objects.get(pk=row["rejected_barcode_id"][0]).name
    #rejection_status = row["rejected_barcode__name"][0]
    status = row["is_pass"][0]
    bin_index = int(row["bin_index"][0])
    sequence_length_sum = int(row["sequence_length"]["sum"])
    read_count = int(row["sequence_length"]["count"])

    uniquekey = "{}_flowcellHistogramSummary_{}_{}_{}_{}_{}".format(flowcell_id, barcode_name,
                                                                    rejection_status, read_type_name, status,bin_index)
    p = redis_instance.pipeline()
    p.multi()
    p.hset(uniquekey, "bin_index", bin_index)
    p.hset(uniquekey, "barcode_name", barcode_name)
    p.hset(uniquekey, "rejection_status", rejection_status)
    p.hset(uniquekey, "read_type_name", read_type_name)
    p.hset(uniquekey, "status", str(status))
    p.hincrby(uniquekey, "read_length", sequence_length_sum)
    p.hincrby(uniquekey, "read_count", read_count)
    p.execute()
    """
    flowcellHistogramSummary, created = FlowcellHistogramSummary.objects.get_or_create(
        flowcell=flowcell,
        barcode_name=barcode_name,
        rejection_status=rejection_status,
        read_type_name=read_type_name,
        status=status,
        bin_index=bin_index,
    )

    flowcellHistogramSummary.read_length += sequence_length_sum
    flowcellHistogramSummary.read_count += read_count

    flowcellHistogramSummary.save()
    """


#s@task(serializer="pickle")
def save_flowcell_summary_barcode_async(row):
    """
    Save flowcell summary barcode rows, applied to a pandas dataframe
    :param flowcell_id: Primary key of the flowcell database entry
    :type flowcell_id: int
    :param row: a pandas dataframe row
    :type row: pandas.core.series.Series
    :return: None
    """
    type_name = FastqReadType.objects.get(pk=row["type_id"][0]).name
    rejection_status = Barcode.objects.get(pk=row["rejected_barcode_id"][0]).name
    flowcell_id = int(row["flowcell_id"][0])
    #flowcell = Flowcell.objects.get(pk=flowcell_id)
    barcode_name = row["barcode_name"][0]
    status = row["is_pass"][0]
    #rejection_status = row["rejected_barcode_name"][0]
    sequence_length_sum = int(row["sequence_length"]["sum"])
    sequence_length_max = int(row["sequence_length"]["max"])
    sequence_length_min = int(row["sequence_length"]["min"])
    quality_average_sum = int(row["quality_average"]["sum"])
    read_count = int(row["sequence_length"]["count"])
    channels = row["channel"]["unique"]
    uniquekey = "{}_flowcellSummaryBarcode_{}_{}_{}_{}".format(flowcell_id, barcode_name,
                                                                    rejection_status, type_name, status)
    prev_channel_list = redis_instance.hget(uniquekey, "channel_presence")
    min_prev = redis_instance.hget(uniquekey, "min_length")
    max_prev = redis_instance.hget(uniquekey, "max_length")
    p = redis_instance.pipeline()
    p.multi()
    p.hset(uniquekey, "barcode_name", barcode_name)
    p.hset(uniquekey, "rejection_status", rejection_status)
    p.hset(uniquekey, "type_name", type_name)
    p.hset(uniquekey, "status", str(status))
    p.hincrby(uniquekey, "total_length", sequence_length_sum)
    p.hincrby(uniquekey, "quality_sum", quality_average_sum)
    p.hincrby(uniquekey, "read_count", read_count)

    if min_prev:
        min_prev = int(min_prev)
        if min_prev > sequence_length_min:
            p.hset(uniquekey, "min_length", sequence_length_min)
    else:
        p.hset(uniquekey, "min_length", sequence_length_min)

    if max_prev:
        max_prev = int(max_prev)
        if max_prev < sequence_length_max:
            p.hset(uniquekey, "max_length", sequence_length_max)
    else:
        p.hset(uniquekey, "max_length", sequence_length_max)

    if prev_channel_list:
        channel_list = list(prev_channel_list)
    else:
        channel_list = list("0" * 3000)

    for c in channels:
        channel_list[int(c) - 1] = "1"

    channel_presence = "".join(channel_list)

    p.hset(uniquekey, "channel_presence", channel_presence)
    p.execute()
    """
    flowcellSummaryBarcode, created = FlowcellSummaryBarcode.objects.get_or_create(
        flowcell=flowcell,
        barcode_name=barcode_name,
        rejection_status=rejection_status,
        read_type_name=type_name,
        status=status,
    )
    flowcellSummaryBarcode.total_length += sequence_length_sum
    flowcellSummaryBarcode.quality_sum += decimal.Decimal(quality_average_sum)
    flowcellSummaryBarcode.read_count += read_count
    if flowcellSummaryBarcode.max_length < sequence_length_max:
        flowcellSummaryBarcode.max_length = sequence_length_max
    if flowcellSummaryBarcode.min_length == 0 or flowcellSummaryBarcode.min_length > sequence_length_min:
        flowcellSummaryBarcode.min_length = sequence_length_min
    channel_list = list(flowcellSummaryBarcode.channel_presence)
    for c in channels:
        channel_list[int(c)  - 1] = "1"
    flowcellSummaryBarcode.channel_presence = "".join(channel_list)
    flowcellSummaryBarcode.channel_count = len(channels)
    flowcellSummaryBarcode.save()
    """


def save_flowcell_summary_barcode(flowcell_id, row):
    """
    Save flowcell summary barcode rows, applied to a pandas dataframe
    :param flowcell_id: Primary key of the flowcell database entry
    :type flowcell_id: int
    :param row: a pandas dataframe row
    :type row: pandas.core.series.Series
    :return: None
    """
    barcode_name = row["barcode__name"][0]
    type_name = row["type__name"][0]
    status = row["is_pass"][0]
    rejection_status = row["rejected_barcode__name"][0]
    sequence_length_sum = row["sequence_length"]["sum"]
    sequence_length_max = row["sequence_length"]["max"]
    sequence_length_min = row["sequence_length"]["min"]
    quality_average_sum = row["quality_average"]["sum"]
    read_count = row["sequence_length"]["count"]
    channels = row["channel"]["unique"]
    try:
        flowcellSummaryBarcode, created = FlowcellSummaryBarcode.objects.get_or_create(
            flowcell_id=flowcell_id,
            barcode_name=barcode_name,
            rejection_status=rejection_status,
            read_type_name=type_name,
            status=status,
        )

        flowcellSummaryBarcode.total_length += sequence_length_sum
        flowcellSummaryBarcode.quality_sum += quality_average_sum
        flowcellSummaryBarcode.read_count += read_count

        if flowcellSummaryBarcode.max_length < sequence_length_max:
            flowcellSummaryBarcode.max_length = sequence_length_max

        if flowcellSummaryBarcode.min_length < sequence_length_min:
            flowcellSummaryBarcode.min_length = sequence_length_min

        channel_list = list(flowcellSummaryBarcode.channel_presence)

        for c in channels:

            channel_list[c - 1] = "1"

        flowcellSummaryBarcode.channel_presence = "".join(channel_list)
        flowcellSummaryBarcode.channel_count = len(channels)
        flowcellSummaryBarcode.save()
    except Exception as e:
        print(e)
        raise

def save_flowcell_statistic_barcode(flowcell_id, row):
    """
    Save flowcell statistic barcode objects into the database row-wise on a pandas dataframe.
    Data for the flowcell statistics on y
    :param flowcell_id: The primary key of the flowcell database entry
    :type flowcell_id: int
    :param row: a pandas dataframe row
    :type row: pandas.core.series.Series
    :return: None
    """

    utc = pytz.utc

    start_time = utc.localize(row["start_time_truncate"][0])
    barcode_name = row["barcode__name"][0]
    rejection_status = row["rejected_barcode__name"][0]
    type_name = row["type__name"][0]
    status = row["is_pass"][0]
    sequence_length_sum = row["sequence_length"]["sum"]
    sequence_length_max = row["sequence_length"]["max"]
    sequence_length_min = row["sequence_length"]["min"]
    quality_average_sum = row["quality_average"]["sum"]
    read_count = row["sequence_length"]["count"]
    channels = row["channel"]["unique"]
    try:
        flowcellStatisticBarcode, created = FlowcellStatisticBarcode.objects.get_or_create(
            flowcell_id=flowcell_id,
            sample_time=start_time,
            barcode_name=barcode_name,
            rejection_status=rejection_status,
            read_type_name=type_name,
            status=status,
        )

        flowcellStatisticBarcode.total_length += sequence_length_sum
        flowcellStatisticBarcode.quality_sum += quality_average_sum
        flowcellStatisticBarcode.read_count += read_count

        if flowcellStatisticBarcode.max_length < sequence_length_max:
            flowcellStatisticBarcode.max_length = sequence_length_max

        if flowcellStatisticBarcode.min_length < sequence_length_min:
            flowcellStatisticBarcode.min_length = sequence_length_min

        channel_list = list(flowcellStatisticBarcode.channel_presence)

        for c in channels:

            channel_list[c - 1] = "1"

            flowcellStatisticBarcode.channel_presence = "".join(channel_list)
            flowcellStatisticBarcode.channel_count = len(channels)

        flowcellStatisticBarcode.save()
    except Exception as e:
        print(e)
        raise

def save_flowcell_histogram_summary(flowcell_id, row):
    """
    Save flowcell histogram barcode objects into the database row-wise on a pandas dataframe.
    Data for the histograms on the chancalc page.
    :param flowcell_id: Primary key of the flowcell database entry
    :type flowcell_id: int
    :param row: a pandas dataframe row
    :type row: pandas.core.series.Series
    :return: None
    """
    barcode_name = row["barcode__name"][0]
    read_type_name = row["type__name"][0]
    rejection_status = row["rejected_barcode__name"][0]
    status = row["is_pass"][0]
    bin_index = row["bin_index"][0]
    sequence_length_sum = row["sequence_length"]["sum"]
    read_count = row["sequence_length"]["count"]
    try:
        flowcellHistogramSummary, created = FlowcellHistogramSummary.objects.get_or_create(
            flowcell_id=flowcell_id,
            barcode_name=barcode_name,
            rejection_status=rejection_status,
            read_type_name=read_type_name,
            status=status,
            bin_index=bin_index,
        )

        flowcellHistogramSummary.read_length += sequence_length_sum
        flowcellHistogramSummary.read_count += read_count

        flowcellHistogramSummary.save()
    except Exception as e:
        print(e)
        raise

def save_flowcell_channel_summary(flowcell_id, row):
    """
    Save flowcell channel summary into the database row-wise on a pandas dataframe. Used for the channel visualisations.
    :param flowcell_id: Primary key of the flowcell database entry
    :type flowcell_id: int
    :param row: a pandas dataframe row
    :type row: pandas.core.series.Series
    :return: None
    """
    channel = row["channel"][0]
    sequence_length_sum = row["sequence_length"]["sum"]
    read_count = row["sequence_length"]["count"]

    flowcellChannelSummary, created = FlowcellChannelSummary.objects.get_or_create(
        flowcell_id=flowcell_id, channel=channel
    )

    flowcellChannelSummary.read_length += sequence_length_sum
    flowcellChannelSummary.read_count += read_count

    flowcellChannelSummary.save()


@receiver(post_save, sender=Flowcell)
def new_flowcell_message(sender, instance=None, created=False, **kwargs):
    """
        If a new Flowcell is created in the database
        Parameters
        ----------
        sender:
            Unknown
        instance: reads.models.Flowcell
            The instance of the flowcell that is being saved.
        created: bool
            Whether the save that is being called is creating the object in the database.
        kwargs: dict
            Key word arguments, if any.

        Returns
        -------
        None

    """
    if created:
        myuser = instance.owner
        new_flowcell_message = Message(
            recipient=myuser,
            sender=myuser,
            flowcell=instance,
            title=f"New flowcell {instance.name} created with sample name {instance.sample_name} "
            f"at {datetime.datetime.now()}.",
        )
        new_flowcell_message.save()


@receiver(post_save, sender=Run)
def new_run_message(sender, instance=None, created=False, **kwargs):
    """
    If a new Run is created in the database
    Parameters
    ----------
    sender:
        Unknown
    instance: reads.models.Run
        The instance of the Run that is being saved.
    created: bool
        Whether the save that is being called is creating the object in the database.
    kwargs: dict
        Key word arguments, if any.

    Returns
    -------
    None

    """
    if created:
        myuser = instance.owner
        new_run_message = Message(
            recipient=myuser,
            sender=myuser,
            title=f"New run {instance.name} created on flowcell {instance.flowcell.name} at {datetime.datetime.now()}.",
            run=instance,
            flowcell=instance.flowcell
        )
        new_run_message.save()


@receiver(post_save, sender=MinionMessage)
def new_minion_message(sender, instance=None, created=False, **kwargs):
    """
    Save a message to the User if the minion message that has come in is a Mux or space warning.

    Parameters
    ----------
    sender:
        Unknown
    instance: reads.models.MinionMessage
        The instance of the MinionMessage that is being saved.
    created: bool
        Whether the save that is being called is creating the object in the database.
    kwargs: dict
        Key word arguments, if any.

    Returns
    -------
    None

    """
    if created:
        myuser = instance.minion.owner
        # Messages sent as Warnings (Severity 2), Messages sent as Errors (Severity 3)
        #print(instance)
        #print(instance.run)
        flowcell = instance.run.flowcell
        queryset = NotificationConditions.objects.filter(
            flowcell=flowcell, completed=False
        )

        if int(instance.severity) > 1:
            # Check if we have notifications switched on for this
            queryset = queryset.filter(notification_type="w/e")
            if queryset:
                title = "{} from computer {} at {}".format(
                    instance.message,
                    instance.minion.computer(),
                    datetime.datetime.now(),
                )
                chunks = wrap(title, 512)
                for chunk in chunks:
                    new_message_message = Message(
                        recipient=myuser,
                        sender=myuser,
                        title=chunk,
                        run=instance.run,
                        flowcell=flowcell
                    )
                    new_message_message.save()
        # Mux messages me thinks
        elif instance.message.startswith("Flow cell"):
            queryset = queryset.filter(notification_type="mux")
            if queryset:
                title = "{} from computer {} at {}".format(
                    instance.message,
                    instance.minion.computer(),
                    datetime.datetime.now(),
                )
                chunks = wrap(title, 512)
                for chunk in chunks:
                    new_message_message = Message(
                        recipient=myuser,
                        sender=myuser,
                        title=chunk,
                        run=instance.run,
                        flowcell=flowcell,
                    )
                    new_message_message.save()
        # Messages sent by minotTour
        elif instance.message.startswith("minoTour"):
            queryset = queryset.filter(notification_type="mino")
            if queryset:
                title = "{} from computer {} at {}".format(
                    instance.message,
                    instance.minion.computer(),
                    datetime.datetime.now(),
                )
                chunks = wrap(title, 512)
                for chunk in chunks:
                    new_message_message = Message(
                        recipient=myuser,
                        sender=myuser,
                        title=chunk,
                        run=instance.run,
                        flowcell=instance.run.flowcell,
                    )
                    new_message_message.save()


### Modifying this to take the results from a bulk save implemented in postgres

@task(serializer="pickle")
def align_reads(fasta_df,flowcell):
    """
    Call minimap2 to align reads to the reference sequence listed on the JobMaster.

    Parameters
    ----------
    fastas: list
        A list of FastqRead objects.
    job_master_id: int
        The PK of the Jobmaster entry.
    fasta_df: pandas.core.frame.DataFrame
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
    #job_master = JobMaster.objects.get(pk=job_master_id)
    job_masters = flowcell.flowcelljobs.filter(job_type__name="Minimap2")

    for job_master in job_masters:
        # Unpack metadata about this job
        reference_info_id = job_master.reference.id

        #
        if not minimap2:
            logger.error("Can not find minimap2 executable - stopping task.")
            return

        # The information about the reference for this alignment
        reference_info = ReferenceInfo.objects.get(pk=reference_info_id)
        # A dictionary with all the results for each chromosome
        chromdict = dict()
        # Get all the chromosomes for this reference
        chromosomes = reference_info.referencelines.all()
        # For each Chromosome in the chromosomes, set the name of the chromosome as key, and the object as the value
        for chromosome in chromosomes:
            chromdict[chromosome.line_name] = chromosome

        # Create the series that contains the read_id and sequence as a correctly formatted fasta string
        fasta_df["fasta"] = ">" + fasta_df["read_id"] + "\n" + fasta_df["sequence"]
        # Create one string formatted correctly for fasta with input data
        fasta_data = "\n".join(list(fasta_df["fasta"]))

        # dict where we will be storing the fastq objects keyed to read ID
        #fastq_dict = {fasta.read_id: fasta for fasta in fastas}
        # Create the minimap command
        # If we have a minimap2 index we want to use it.
        if reference_info.minimap2_index_file_location:
            referencefile = reference_info.minimap2_index_file_location
        else:
            referencefile = reference_info.filename
        cmd = '{} -x map-ont -t 4 --secondary=no {} -'.format(
            minimap2,
            os.path.join(reference_location, referencefile)
        )

        logger.info("Flowcell id: {} - Calling minimap - {}".format(flowcell.id, cmd))
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
            logger.info(err)

        logger.info("Flowcell id: {} - Finished minimap - {}".format(flowcell.id, cmd))
        # Get the results of the Paf output
        paf = out
        # Split the lines of the paf file
        pafdata = paf.splitlines()

        logger.info(
            "Flowcell id: {} - Found {} paf records".format(flowcell.id, len(pafdata))
        )
        # If we have paf data
        if pafdata:

            from django.db import close_old_connections
            close_old_connections()

            incstore=dict()

            bulk_paf = []

            bulk_paf_rough = []
            # For each line in the paf data
            for line in pafdata:
                # Strip the new line char from the file
                line = line.strip("\n")
                # Split the line on the tab character
                record = line.split("\t")
                # Get the fastqread object to use in the Foreign key field below
                #fastq_read = fastq_dict[record[0]]
                fastq_read = fasta_df[fasta_df["read_id"]==record[0]]
                #newpaf = PafStore(job_master=job_master,read=FastqRead.objects.get(pk=fastq_read["id"].values[0]))
                newpaf = PafStore(job_master=job_master,read=FastqRead.objects.get(read_id=fastq_read["read_id"].values[0],flowcell=flowcell))
                #fastq_read.barcode_name = "All reads" if fastq_read.barcode_name == "No barcode" else fastq_read.barcode_name
                barcode_name = fastq_read["barcode_name"].values[0]
                barcode_name = "All reads" if barcode_name == "No barcode" else barcode_name
                fastq_read["barcode_name"]=barcode_name

                newpaf.reference = reference_info
                chromosome_line_name = chromdict[record[5]]
                # Set the fields on the PafStore
                newpaf.qsn = (
                    #fastq_read.read_id
                    record[0]
                )  # models.CharField(max_length=256)#1	string	Query sequence name
                newpaf.qsl = int(
                    record[1]
                )  # models.IntegerField()#2	int	Query sequence length
                newpaf.qs = int(
                    record[2]
                )  # models.IntegerField()#3	int	Query start (0-based)
                newpaf.qe = int(
                    record[3]
                )  # models.IntegerField()#4	int	Query end (0-based)
                newpaf.rs = record[
                    4
                ]  # models.CharField(max_length=1)#5	char	Relative strand: "+" or "-"
                newpaf.tsn = chromosome_line_name  # models.CharField(max_length=256)#6	string	Target sequence name
                newpaf.tsl = int(
                    record[6]
                )  # models.IntegerField()#7	int	Target sequence length
                newpaf.ts = int(
                    record[7]
                )  # models.IntegerField()#8	int	Target start on original strand (0-based)
                newpaf.te = int(
                    record[8]
                )  # models.IntegerField()#9	int	Target end on original strand (0-based)
                newpaf.nrm = int(
                    record[9]
                )  # models.IntegerField()#10	int	Number of residue matches
                newpaf.abl = int(
                    record[10]
                )  # models.IntegerField()#11	int	Alignment block length
                newpaf.mq = int(
                    record[11]
                )  # models.IntegerField()#12	int	Mapping quality (0-255; 255 for missing)
                # Set the fields on the paf rough cov object


                if fastq_read['type_id'].values[0] not in incstore.keys():
                    incstore[fastq_read['type_id'].values[0]]=dict()

                if fastq_read["barcode_name"].values[0] not in incstore[fastq_read['type_id'].values[0]].keys():
                    incstore[fastq_read['type_id'].values[0]][fastq_read["barcode_name"].values[0]]=dict()

                if fastq_read['is_pass'].values[0] not in incstore[fastq_read['type_id'].values[0]][fastq_read["barcode_name"].values[0]].keys():
                    incstore[fastq_read['type_id'].values[0]][fastq_read["barcode_name"].values[0]][fastq_read['is_pass'].values[0]]=dict()

                if reference_info not in incstore[fastq_read['type_id'].values[0]][fastq_read["barcode_name"].values[0]][fastq_read['is_pass'].values[0]].keys():
                    incstore[fastq_read['type_id'].values[0]][fastq_read["barcode_name"].values[0]][fastq_read['is_pass'].values[0]][reference_info]=dict()

                if chromosome_line_name not in incstore[fastq_read['type_id'].values[0]][fastq_read["barcode_name"].values[0]][fastq_read['is_pass'].values[0]][reference_info].keys():
                    incstore[fastq_read['type_id'].values[0]][fastq_read["barcode_name"].values[0]][fastq_read['is_pass'].values[0]][reference_info][chromosome_line_name]=dict()
                    incstore[fastq_read['type_id'].values[0]][fastq_read["barcode_name"].values[0]][fastq_read['is_pass'].values[0]][reference_info][chromosome_line_name]['start']=defaultdict(int)
                    incstore[fastq_read['type_id'].values[0]][fastq_read["barcode_name"].values[0]][fastq_read['is_pass'].values[0]][reference_info][chromosome_line_name]['end']=defaultdict(int)

                incstore[fastq_read['type_id'].values[0]][fastq_read["barcode_name"].values[0]][fastq_read['is_pass'].values[0]][reference_info][chromosome_line_name]['start'][int(record[7])]+=1
                incstore[fastq_read['type_id'].values[0]][fastq_read["barcode_name"].values[0]][fastq_read['is_pass'].values[0]][reference_info][chromosome_line_name]['end'][int(record[8])]-=1


                # newpafstart = PafRoughCov(
                #     job_master=job_master,
                #     flowcell=flowcell,
                #     read_type_id=fastq_read['type_id'].values[0],
                #     barcode_name=fastq_read["barcode_name"].values[0],
                #     is_pass=fastq_read['is_pass'].values[0],
                # )
                # newpafend = PafRoughCov(
                #     job_master=job_master,
                #     flowcell=flowcell,
                #     read_type_id=fastq_read['type_id'].values[0],
                #     barcode_name=fastq_read["barcode_name"].values[0],
                #     is_pass=fastq_read['is_pass'].values[0],
                # )

                #newpafstart.reference = reference_info
                #newpafstart.chromosome = chromosome_line_name
                #newpafstart.p = int(record[7])
                #newpafstart.i = 1
                # Set the end fields
                #newpafend.reference = reference_info
                #newpafend.chromosome = chromosome_line_name
                #newpafend.p = int(record[8])
                #newpafend.i = -1
                #bulk_paf_rough.append(newpafstart)
                #bulk_paf_rough.append(newpafend)
                bulk_paf.append(newpaf)

            for read_type_id in incstore:
                for barcode_name in incstore[read_type_id]:
                    for is_pass in incstore[read_type_id][barcode_name]:
                        for reference_info in incstore[read_type_id][barcode_name][is_pass]:
                            for chromosome_line_name in incstore[read_type_id][barcode_name][is_pass][reference_info]:
                                for p in incstore[read_type_id][barcode_name][is_pass][reference_info][chromosome_line_name]['start']:
                                    newpafstart = PafRoughCov(
                                        job_master=job_master,
                                        flowcell=flowcell,
                                        read_type_id=read_type_id,
                                        barcode_name=barcode_name,
                                        is_pass=is_pass,
                                    )
                                    newpafstart.reference = reference_info
                                    newpafstart.chromosome = chromosome_line_name
                                    newpafstart.p = p
                                    newpafstart.i = incstore[read_type_id][barcode_name][is_pass][reference_info][chromosome_line_name]['start'][p]
                                    bulk_paf_rough.append(newpafstart)

                                for p in incstore[read_type_id][barcode_name][is_pass][reference_info][chromosome_line_name]['end']:
                                    newpafend = PafRoughCov(
                                        job_master=job_master,
                                        flowcell=flowcell,
                                        read_type_id=read_type_id,
                                        barcode_name=barcode_name,
                                        is_pass=is_pass,
                                    )
                                    newpafend.reference = reference_info
                                    newpafend.chromosome = chromosome_line_name
                                    newpafend.p = p
                                    newpafend.i = \
                                    incstore[read_type_id][barcode_name][is_pass][reference_info][chromosome_line_name][
                                        'end'][p]
                                    bulk_paf_rough.append(newpafend)


            PafStore.objects.bulk_create(bulk_paf, batch_size=1000)
            PafRoughCov.objects.bulk_create(bulk_paf_rough, batch_size=1000)

            paf_store_list = PafStore.objects.filter(job_master=job_master).values(
                "read__barcode_name",
                "tsn__line_name",
                "tsn__id",
                "qsn",
                "qs",
                "qe",
                "read__type__name",
            )
            # Create a DataFrame of the PafStore results from above, plus previous
            paf_store_df = pd.DataFrame.from_records(paf_store_list)

            paf_store_df["read__barcode_name"][
                paf_store_df["read__barcode_name"] == "No barcode"
            ] = "All reads"

            if np.unique(paf_store_df["read__barcode_name"].values).size > 1:

                paf_store_df_all_reads = paf_store_df.copy(deep=True)

                paf_store_df_all_reads["read__barcode_name"] = "All reads"

                paf_store_df_all_reads["length"] = (
                    paf_store_df_all_reads["qe"] - paf_store_df_all_reads["qs"]
                )

                paf_store_df_all_reads_gb = paf_store_df_all_reads.groupby(
                    ["read__barcode_name", "tsn__line_name", "read__type__name"]
                ).agg({"qsn": ["unique"], "length": ["sum"]})
                # Apply a saving function down the dataframe to save objects
                paf_store_df_all_reads_gb.reset_index().apply(
                    lambda row: save_paf_store_summary(job_master, row), axis=1
                )

            # Minus the query start from the query end
            paf_store_df["length"] = paf_store_df["qe"] - paf_store_df["qs"]
            # Do some maths?
            paf_store_gb = paf_store_df.groupby(
                ["read__barcode_name", "tsn__line_name", "read__type__name"]
            ).agg({"qsn": ["unique"], "length": ["sum"]})
            # Apply a saving function down the dataframe to save objects
            paf_store_gb.reset_index().apply(
                lambda row: save_paf_store_summary(job_master, row), axis=1
            )



def save_paf_store_summary(job_master, row):
    """
    Save the Paf Store Summaries
    :param job_master: The job master object
    :param row: The row of the dataframe
    :return:
    """
    barcode_name = row["read__barcode_name"][0]
    reference_line_name = row["tsn__line_name"][0]
    total_length = row["length"]["sum"]
    read_list = row["qsn"]["unique"]
    read_type = row["read__type__name"][0]
    # The reference line for this object
    reference_line = ReferenceLine.objects.filter(
        reference=job_master.reference
    ).filter(line_name=reference_line_name)[0]
    read_type = FastqReadType.objects.get(name=read_type)
    paf_summary_cov, created = PafSummaryCov.objects.get_or_create(
        job_master_id=job_master.id,
        barcode_name=barcode_name,
        chromosome=reference_line,
        reference_line_length=reference_line.chromosome_length,
        read_type=read_type
    )
    paf_summary_cov.total_length = total_length
    paf_summary_cov.read_count = len(read_list)
    paf_summary_cov.save()
