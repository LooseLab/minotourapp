import pytz

from reads.models import FlowcellStatisticBarcode, FlowcellSummaryBarcode, FlowcellHistogramSummary, FlowcellChannelSummary


def save_flowcell_summary_barcode(flowcell_id, row):

    barcode_name = row['barcode__name'][0]
    type_name = row['type__name'][0]
    status = row['is_pass'][0]
    sequence_length_sum = row['sequence_length']['sum']
    sequence_length_max = row['sequence_length']['max']
    sequence_length_min = row['sequence_length']['min']
    quality_average_sum = row['quality_average']['sum']
    read_count = row['sequence_length']['count']
    channels = row['channel']['unique']

    flowcellSummaryBarcode, created = FlowcellSummaryBarcode.objects.get_or_create(
        flowcell_id=flowcell_id,
        barcode_name=barcode_name,
        read_type_name=type_name,
        status=status)

    flowcellSummaryBarcode.total_length += sequence_length_sum
    flowcellSummaryBarcode.quality_sum += quality_average_sum
    flowcellSummaryBarcode.read_count += read_count

    if flowcellSummaryBarcode.max_length < sequence_length_max:
        flowcellSummaryBarcode.max_length = sequence_length_max

    if flowcellSummaryBarcode.min_length < sequence_length_min:
        flowcellSummaryBarcode.min_length = sequence_length_min

    channel_list = list(flowcellSummaryBarcode.channel_presence)

    for c in channels:

        channel_list[c - 1] = '1'

    flowcellSummaryBarcode.channel_presence = ''.join(channel_list)
    flowcellSummaryBarcode.channel_count = len(channels)
    flowcellSummaryBarcode.save()


def save_flowcell_statistic_barcode(flowcell_id, row):

    utc = pytz.utc

    start_time = utc.localize(row['start_time_truncate'][0])
    barcode_name = row['barcode__name'][0]
    type_name = row['type__name'][0]
    status = row['is_pass'][0]
    sequence_length_sum = row['sequence_length']['sum']
    sequence_length_max = row['sequence_length']['max']
    sequence_length_min = row['sequence_length']['min']
    quality_average_sum = row['quality_average']['sum']
    read_count = row['sequence_length']['count']
    channels = row['channel']['unique']

    flowcellStatisticBarcode, created = FlowcellStatisticBarcode.objects.get_or_create(
        flowcell_id=flowcell_id,
        sample_time=start_time,
        barcode_name=barcode_name,
        read_type_name=type_name,
        status=status)

    flowcellStatisticBarcode.total_length += sequence_length_sum
    flowcellStatisticBarcode.quality_sum += quality_average_sum
    flowcellStatisticBarcode.read_count += read_count

    if flowcellStatisticBarcode.max_length < sequence_length_max:
        flowcellStatisticBarcode.max_length = sequence_length_max

    if flowcellStatisticBarcode.min_length < sequence_length_min:
        flowcellStatisticBarcode.min_length = sequence_length_min

    channel_list = list(flowcellStatisticBarcode.channel_presence)

    for c in channels:

        channel_list[c - 1] = '1'

        flowcellStatisticBarcode.channel_presence = ''.join(channel_list)
        flowcellStatisticBarcode.channel_count = len(channels)

    flowcellStatisticBarcode.save()


def save_flowcell_histogram_summary(flowcell_id, row):

    barcode_name = row['barcode__name'][0]
    read_type_name = row['type__name'][0]
    status = row['is_pass'][0]
    bin_index = row['bin_index'][0]
    sequence_length_sum = row['sequence_length']['sum']
    read_count = row['sequence_length']['count']

    flowcellHistogramSummary, created = FlowcellHistogramSummary.objects.get_or_create(
        flowcell_id=flowcell_id,
        barcode_name=barcode_name,
        read_type_name=read_type_name,
        status=status,
        bin_index=bin_index)

    flowcellHistogramSummary.read_length += sequence_length_sum
    flowcellHistogramSummary.read_count += read_count

    flowcellHistogramSummary.save()


def save_flowcell_channel_summary(flowcell_id, row):

    channel = row['channel'][0]
    sequence_length_sum = row['sequence_length']['sum']
    read_count = row['sequence_length']['count']

    flowcellChannelSummary, created = FlowcellChannelSummary.objects.get_or_create(
        flowcell_id=flowcell_id,
        channel=channel)

    flowcellChannelSummary.read_length += sequence_length_sum
    flowcellChannelSummary.read_count += read_count

    flowcellChannelSummary.save()
