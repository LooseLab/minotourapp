from reads.models import RunSummaryBarcode, RunStatisticBarcode, HistogramSummary, ChannelSummary, FlowcellStatisticBarcode, FlowcellSummaryBarcode, FlowcellHistogramSummary, FlowcellChannelSummary


def save_runsummarybarcode(run_id, row):

    barcode_id = row['barcode_id'][0]
    type_id = row['type_id'][0]
    is_pass = row['is_pass'][0]
    sequence_length_sum = row['sequence_length']['sum']
    sequence_length_max = row['sequence_length']['max']
    sequence_length_min = row['sequence_length']['min']
    quality_average_sum = row['quality_average']['sum']
    read_count = row['sequence_length']['count']
    channels = row['channel']['unique']

    runSummaryBarcode, created = RunSummaryBarcode.objects.get_or_create(
        run_id = run_id,
        barcode_id = barcode_id,
        type_id = type_id,
        is_pass = is_pass)

    runSummaryBarcode.total_length += sequence_length_sum
    runSummaryBarcode.quality_sum += quality_average_sum
    runSummaryBarcode.read_count += read_count

    if runSummaryBarcode.max_length < sequence_length_max:
        runSummaryBarcode.max_length = sequence_length_max

    if runSummaryBarcode.min_length < sequence_length_min:
        runSummaryBarcode.min_length = sequence_length_min

    channel_list = list(runSummaryBarcode.channel_presence)

    for c in channels:

        channel_list[c - 1] = '1'

    runSummaryBarcode.channel_presence = ''.join(channel_list)
    runSummaryBarcode.channel_count = len(channels)
    runSummaryBarcode.save()


def save_flowcellsummarybarcode(flowcell_id, row):

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


def save_runstatisticbarcode(run_id, row):

    start_time = row['start_time'][0]
    barcode_id = row['barcode_id'][0]
    type_id = row['type_id'][0]
    is_pass = row['is_pass'][0]
    sequence_length_sum = row['sequence_length']['sum']
    sequence_length_max = row['sequence_length']['max']
    sequence_length_min = row['sequence_length']['min']
    quality_average_sum = row['quality_average']['sum']
    read_count = row['sequence_length']['count']
    channels = row['channel']['unique']

    runStatisticBarcode, created = RunStatisticBarcode.objects.get_or_create(
        run_id = run_id,
        sample_time = start_time,
        barcode_id = barcode_id,
        type_id = type_id,
        is_pass = is_pass)

    runStatisticBarcode.total_length += sequence_length_sum
    runStatisticBarcode.quality_sum += quality_average_sum
    runStatisticBarcode.read_count += read_count

    if runStatisticBarcode.max_length < sequence_length_max:
        runStatisticBarcode.max_length = sequence_length_max

    if runStatisticBarcode.min_length < sequence_length_min:
        runStatisticBarcode.min_length = sequence_length_min

    channel_list = list(runStatisticBarcode.channel_presence)

    for c in channels:

        channel_list[c - 1] = '1'

        runStatisticBarcode.channel_presence = ''.join(channel_list)
        runStatisticBarcode.channel_count = len(channels)

    runStatisticBarcode.save()


def save_flowcellstatisticbarcode(flowcell_id, row):

    start_time = row['start_time_truncate'][0]
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


def save_histogramsummary(run_id, row):

    barcode_id = row['barcode_id'][0]
    read_type_id = row['type_id'][0]
    is_pass = row['is_pass'][0]
    bin_index = row['bin_index'][0]
    sequence_length_sum = row['sequence_length']['sum']
    read_count = row['sequence_length']['count']

    runHistogramSummary, created = HistogramSummary.objects.get_or_create(
        run_id = run_id,
        barcode_id = barcode_id,
        read_type_id = read_type_id,
        is_pass = is_pass,
        bin_index = bin_index)

    runHistogramSummary.read_length += sequence_length_sum
    runHistogramSummary.read_count += read_count

    runHistogramSummary.save()


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


def save_channelsummary(run_id, row):

    channel = row['channel'][0]
    sequence_length_sum = row['sequence_length']['sum']
    read_count = row['sequence_length']['count']

    runChannelSummary, created = ChannelSummary.objects.get_or_create(
        run_id = run_id,
        channel = channel)

    runChannelSummary.read_length += sequence_length_sum
    runChannelSummary.read_count += read_count

    runChannelSummary.save()


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

