"""
Services to help run some of the tasks in the web app. Why is this not stored there?
"""
import pytz
import datetime
import pandas as pd
import numpy as np
import decimal

from celery.task import task
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.db import IntegrityError, transaction

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

from communication.models import Message, NotificationConditions

import redis

redis_instance = redis.StrictRedis(host="127.0.0.1",
                                  port=6379, db=0,decode_responses=True)


# @task(rate_limit="2/m")

@task()
def save_reads_bulk(reads):
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
    try:
        FastqRead.objects.bulk_create(reads_list)
        update_flowcell.delay(reads_list)
    except Exception as e:
        print(e)
        return str(e)




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
    readDF = pd.DataFrame([o.__dict__ for o in reads_list]).drop(['_state'], axis=1)
    readDF['channel'] = readDF['channel'].astype(int)
    readDF['sequence_length'] = readDF['sequence_length'].astype(int)

    results = readDF.groupby(['flowcell_id']).agg({"read_id": "count", "channel": "max","sequence_length":"sum"})
    results.reset_index().apply(
        lambda row: update_flowcell_counts(row), axis=1
    )


    ## We can assume that we will always have only a single flowcell in a batch of reads to process due to the way that minFQ behaves.
    ## Todo: catch problems with multiple flowcells in a single transaction

    readDF["status"] = np.where(
        readDF["is_pass"] == False, "Fail", "Pass"
    )

    readDF["start_time_truncate"] = np.array(
        readDF["start_time"], dtype="datetime64[m]"
    )

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
    redis_instance.set("{}_has_fastq".format(int(row["flowcell_id"])),1)
    redis_instance.incrby("{}_number_reads".format(int(row["flowcell_id"])),int(row["read_id"]))
    redis_instance.incrby("{}_total_read_length".format(int(row["flowcell_id"])),int(row["sequence_length"]))
    saved_max_channel = redis_instance.get("{}_max_channel".format(int(row["flowcell_id"])))
    if saved_max_channel:
            if (row["channel"])>int(saved_max_channel):
                redis_instance.set("{}_max_channel".format(int(row["flowcell_id"])),int(row["channel"]))
    else:
        redis_instance.set("{}_max_channel".format(int(row["flowcell_id"])),int(row["channel"]))
    """
    flowcell = (Flowcell.objects.select_for_update().filter(pk=int(row["flowcell_id"]))[0])
    # We want to try with this being atomic to make sure only one task works on the row at once.
    # Get the flowcell we are working on:
    #flowcell = Flowcell.objects.get(pk=row["flowcell_id"])
    if not flowcell.has_fastq:
        flowcell.has_fastq = True
    #print (row)
    #this is actually the count of the number of read_ids
    flowcell.number_reads += row["read_id"]
    #update the total read length on the flowcell
    flowcell.total_read_length += row["sequence_length"]
    #flowcell.save()
    flowcell.average_read_length = flowcell.total_read_length/flowcell.number_reads
    if int(row['channel']) > flowcell.max_channel:
        flowcell.max_channel = int(row['channel'])
        if row['channel']:
            if int(row['channel']) > 512:
                flowcell.size = 3000
            elif int(row['channel']) > 126:
                flowcell.size = 512
            else:
                flowcell.size = 126
        else:
            flowcell.size = 512
    flowcell.save()
    return
    """


@task()
def update_chan_calc(readDF):
    pass

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
    p.hincrby(uniquekey,"total_length",sequence_length_sum)
    p.hincrby(uniquekey,"quality_sum",quality_average_sum)
    p.hincrby(uniquekey,"read_count",read_count)
    if min_prev:
        min_prev=int(min_prev)
        if min_prev > sequence_length_min:
            p.hset(uniquekey,"min_length",sequence_length_min)
    else:
        p.hset(uniquekey, "min_length", sequence_length_min)
    if max_prev:
        max_prev=int(max_prev)
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
    flowcell = Flowcell.objects.get(pk=flowcell_id)
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
