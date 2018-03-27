from .models import RunSummaryBarcode, RunStatisticBarcode, HistogramSummary, ChannelSummary


def update_flowcelldetails(flowcell_id):

    flowcell = Flowcell.objects.get(pk=flowcell_id)

    flowcell_details = flowcell.flowcell_details

    # check if flowcell has barcoded run

    flowcellruns = flowcell.flowcelldetails

    flowcell_details.has_barcode = False

    for flowcellrun in flowcellruns.all():
        run = flowcellrun.run
        if run.is_barcoded:
            flowcell_details.has_barcode = True
            break

    flowcell_details.save()


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

    """print("barcode_id: {}".format(barcode_id))
    print("type_id: {}".format(type_id))
    print("is_pass: {}".format(is_pass))
    print("sequence_length min: {}".format(row['sequence_length']['min']))
    print("sequence_length max: {}".format(row['sequence_length']['max']))
    print("sequence_length sum: {}".format(row['sequence_length']['sum']))
    print("count: {}".format(row['sequence_length']['count']))
    print("quality average sum: {}".format(row['quality_average']['sum']))
    print("channels: ")
    print(channels)
    print("channel list: {}".format(''.join(channel_list)))"""

    runSummaryBarcode, created = RunSummaryBarcode.objects.get_or_create(
        run_id_id = run_id,
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
        run_id_id = run_id,
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
        run_id_id = run_id,
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


