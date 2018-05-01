from reads.models import RunSummaryBarcode, RunStatisticBarcode, HistogramSummary, ChannelSummary, Run, GroupBarcode


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
    runSummaryBarcode.save()


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

    runStatisticBarcode.save()


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

