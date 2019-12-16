"""
Services to help run some of the tasks in the web app. Why is this not stored there?
"""
import pytz
import datetime
from django.db.models.signals import post_save
from django.dispatch import receiver

from reads.models import (
    FlowcellStatisticBarcode,
    FlowcellSummaryBarcode,
    FlowcellHistogramSummary,
    FlowcellChannelSummary,
    Flowcell, JobMaster, JobType, Run)

from communication.models import Message


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
    Log a message in the system if a new flowcell is created.
    :param sender:
    :param instance:
    :param created:
    :param kwargs:
    :return:
    """

    if created:
        myuser = instance.owner
        new_flowcell_message = Message(recipient=myuser,sender=myuser,title="New flowcell {} created with sample name {} at {}.".format(instance.name, instance.sample_name, datetime.datetime.now()))
        new_flowcell_message.save()


@receiver(post_save, sender=Run)
def new_run_message(sender, instance=None, created=False, **kwargs):
    """

    :param sender:
    :param instance:
    :param created:
    :param kwargs:
    :return:
    """
    if created:
        myuser = instance.owner
        new_run_message = Message(recipient=myuser,sender=myuser,title="New run {} created on flowcell {} at {}.".format(instance.name,instance.flowcell.name,datetime.datetime.now()))
        new_run_message.save()

# @receiver(post_save, sender=Flowcell)
# def create_flowcell_jobs(sender, instance=None, created=False, **kwargs):
#     """
#     Create the flowcell jobs automatically on flowcell creation
#     Parameters
#     ----------
#     sender
#     instance: reads.models.Flowcell
#         The flowcell database model
#     created: bool
#         Whether the flowcell has been created
#     kwargs
#
#     Returns
#     -------
#
#     """
#
#     if created:
#
#         job_name_list = ['ChanCalc', 'UpdateFlowcellDetails']
#
#         for job_name in job_name_list:
#
#             JobMaster.objects.get_or_create(
#                 job_type=JobType.objects.get(name=job_name),
#                 flowcell=instance,
#             )
