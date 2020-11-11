"""
Celery task to update flowcell metadata
"""

import numpy as np
from celery import chain
from celery.utils.log import get_task_logger
from dateutil import parser

from minotourapp.celery import app
from minotourapp.redis import redis_instance
from reads.models import (
    FlowcellHistogramSummary,
    FlowcellStatisticBarcode,
    FlowcellSummaryBarcode,
    JobMaster,
    MinionRunInfo,
    Run,
    FlowcellChannelSummary, Barcode,
)
from reads.tasks.redis_tasks_functions import (
    get_values_and_delete_key_redis_hash,
    get_values_and_delete_redis_key,
)

logger = get_task_logger(__name__)


@app.task
def calculate_flowcell_summary_barcode(flowcell_id):
    """
    Flowcell summary barcode, get all keys related to the flowcell summary barcode values of this flowcell
    from redis and update them
    Parameters
    ----------
    flowcell_id: int

    Returns
    -------

    """
    # flowcell = Flowcell.objects.get(pk=flowcell_id)
    keys = redis_instance.scan_iter("{}_flowcellSummaryBarcode_*".format(flowcell_id))
    for key in keys:
        result = get_values_and_delete_key_redis_hash(redis_instance, key)
        barcode_name = result["barcode_name"]
        rejection_status = result["rejection_status"]
        type_name = result["type_name"]
        status = result["status"]
        sequence_length_sum = int(result["total_length"])
        quality_average_sum = int(result["quality_sum"])
        read_count = int(result["read_count"])
        sequence_length_max = int(result["max_length"])
        sequence_length_min = int(result["min_length"])
        new_channel_list = list(result["channel_presence"])

        (
            flowcell_summary_barcode,
            created,
        ) = FlowcellSummaryBarcode.objects.get_or_create(
            flowcell_id=flowcell_id,
            barcode_name=barcode_name,
            rejection_status=rejection_status,
            read_type_name=type_name,
            status=status,
        )

        flowcell_summary_barcode.total_length += sequence_length_sum
        flowcell_summary_barcode.quality_sum += quality_average_sum
        flowcell_summary_barcode.read_count += read_count
        if flowcell_summary_barcode.max_length < sequence_length_max:
            flowcell_summary_barcode.max_length = sequence_length_max
        if (
            flowcell_summary_barcode.min_length == 0
            or flowcell_summary_barcode.min_length > sequence_length_min
        ):
            flowcell_summary_barcode.min_length = sequence_length_min

        channel_list = list(flowcell_summary_barcode.channel_presence)
        fusedlist = np.bitwise_or(
            np.array(channel_list, dtype=int), np.array(new_channel_list, dtype=int)
        )

        flowcell_summary_barcode.channel_presence = ("").join(
            list(fusedlist.astype(str))
        )
        flowcell_summary_barcode.channel_count = np.count_nonzero(fusedlist == 1)

        flowcell_summary_barcode.save()
    return None


@app.task
def calculate_flowcell_histogram_summary(flowcell_id):
    """
    Flowcell histogram summary. Get all keys related to the flowcell histogram values of this flowcell
    from redis and update them
    Parameters
    ----------
    flowcell_id: int

    Returns
    -------

    """
    # flowcell = Flowcell.objects.get(pk=flowcell_id)
    keys = redis_instance.scan_iter("{}_flowcellHistogramSummary_*".format(flowcell_id))
    for key in keys:
        result = get_values_and_delete_key_redis_hash(redis_instance, key)
        bin_index = int(result["bin_index"])
        barcode_name = result["barcode_name"]
        rejection_status = result["rejection_status"]
        read_type_name = result["read_type_name"]
        status = result["status"]
        sequence_length_sum = int(result["read_length"])
        read_count = int(result["read_count"])

        (
            flowcell_histogram_summary,
            created,
        ) = FlowcellHistogramSummary.objects.get_or_create(
            flowcell_id=flowcell_id,
            barcode_name=barcode_name,
            rejection_status=rejection_status,
            read_type_name=read_type_name,
            status=status,
            bin_index=bin_index,
        )

        flowcell_histogram_summary.read_length += sequence_length_sum
        flowcell_histogram_summary.read_count += read_count

        flowcell_histogram_summary.save()

    return flowcell_id


@app.task
def calculate_flowcell_statistic_barcodes(flowcell_id):
    """
    Flowcell statistic Barcode. Get all keys related to the flowcell summary barcode values of this flowcell
    from redis and update them
    Parameters
    ----------
    flowcell_id: int

    Returns
    -------

    """
    # flowcell = Flowcell.objects.get(pk=flowcell_id)
    keys = redis_instance.scan_iter("{}_flowcellStatisticBarcode_*".format(flowcell_id))
    for key in keys:
        result = get_values_and_delete_key_redis_hash(redis_instance, key)
        start_time = parser.parse(result["sample_time"])
        barcode_name = result["barcode_name"]
        rejection_status = result["rejection_status"]
        type_name = result["type_name"]
        status = result["status"]
        sequence_length_sum = int(result["total_length"])
        quality_average_sum = int(result["quality_sum"])
        read_count = int(result["read_count"])
        sequence_length_max = int(result["max_length"])
        sequence_length_min = int(result["min_length"])
        new_channel_list = list(result["channel_presence"])

        (
            flowcellStatisticBarcode,
            created,
        ) = FlowcellStatisticBarcode.objects.get_or_create(
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
        fusedlist = np.bitwise_or(
            np.array(channel_list, dtype=int), np.array(new_channel_list, dtype=int)
        )

        flowcellStatisticBarcode.channel_presence = ("").join(
            list(fusedlist.astype(str))
        )
        flowcellStatisticBarcode.channel_count = np.count_nonzero(fusedlist == 1)

        flowcellStatisticBarcode.save()

    return flowcell_id


@app.task
def update_flowcell_details(job_master_id):
    """
    Update flowcell details by pulling Using the MinionRunInfo if it exists. If base-called data present,
     pull down aggregated numbers from Redis and create database emtries for them.
    Parameters
    ----------
    job_master_id: int
        The primary key of the job master record for this task.
    Returns
    -------
    None
    """

    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = True
    job_master.save()
    flowcell = job_master.flowcell

    logger.info(
        f"Flowcell id: {flowcell.id} - Updating details of flowcell {flowcell.name}"
    )

    #
    # Get the first MinIONRunStatus for a particular flowcell - but we want to ignore platform QCs
    #
    minion_run_status_first = (
        MinionRunInfo.objects.filter(run_id__flowcell=flowcell)
        .exclude(experiment_type="platform_qc")
        .order_by("minKNOW_start_time")
        .first()
    )
    #
    # If the MinIONRunStatus exists, than update start time
    #
    if minion_run_status_first:
        flowcell.start_time = minion_run_status_first.minKNOW_start_time

    # Get minionRunInfos looking for sample name
    minion_run_info_list = MinionRunInfo.objects.filter(
        run_id__flowcell=flowcell
    ).exclude(experiment_type="platform_qc")

    if minion_run_info_list.count() > 0:
        for minion_run_status in minion_run_info_list:
            if minion_run_status.minKNOW_sample_name != "undefined":
                flowcell.sample_name = minion_run_status_first.minKNOW_sample_name
                break
    # No live data? No problem!
    else:
        for run in Run.objects.filter(flowcell=flowcell).exclude(name="mux scan"):
            if run.name != "undefined":
                flowcell.sample_name = run.name
                if run.start_time is not None:
                    if run.start_time <= flowcell.start_time:
                        flowcell.start_time = run.start_time
                break

    # Get flowcell max_channel number.
    has_fastq = get_values_and_delete_redis_key(
        redis_instance, f"{flowcell.id}_has_fastq"
    )
    if has_fastq:
        flowcell.has_fastq = has_fastq
        total_read_length = get_values_and_delete_redis_key(
            redis_instance, f"{flowcell.id}_total_read_length"
        )
        if total_read_length:
            flowcell.total_read_length += int(total_read_length)

        number_reads = get_values_and_delete_redis_key(
            redis_instance, f"{flowcell.id}_number_reads"
        )
        if number_reads:
            flowcell.number_reads += int(number_reads)

        if number_reads and total_read_length:
            flowcell.average_read_length = int(total_read_length) / int(number_reads)

        max_channel = get_values_and_delete_redis_key(
            redis_instance, f"{flowcell.id}_max_channel"
        )
        if max_channel:
            max_channel = int(max_channel)
            if max_channel > flowcell.max_channel:
                flowcell.max_channel = max_channel
                if max_channel > 512:
                    flowcell.size = 3000
                elif max_channel > 126:
                    flowcell.size = 512
                else:
                    flowcell.size = 126
        else:
            flowcell.size = 512

        # Now to try and update all the channel values that we need to update...

        for channel in range(1, flowcell.max_channel + 1):
            read_count = get_values_and_delete_redis_key(
                redis_instance, f"{flowcell.id}_{channel}_read_count"
            )
            pore_yield = get_values_and_delete_redis_key(
                redis_instance, f"{flowcell.id}_{channel}_pore_yield"
            )
            if read_count or pore_yield:
                (
                    flowcellChannelSummary,
                    created,
                ) = FlowcellChannelSummary.objects.get_or_create(
                    flowcell=flowcell, channel=channel
                )
                if read_count:
                    flowcellChannelSummary.read_count += int(read_count)
                if pore_yield:
                    flowcellChannelSummary.read_length += int(pore_yield)
                flowcellChannelSummary.save()

        #Now try and update all the flowcellsummarystatistics

        res = chain(
            calculate_flowcell_statistic_barcodes.s(flowcell.id),
            calculate_flowcell_histogram_summary.s(),
            calculate_flowcell_summary_barcode.s(),
        ).apply_async()

    flowcell.number_runs = Run.objects.filter(flowcell=flowcell).count()

    barcode_count = (
        Barcode.objects.filter(run_id__flowcell=flowcell).values("name").distinct().count() - 1
    )
    flowcell.number_barcodes = barcode_count
    flowcell.save()
    logger.info(
        "Flowcell id: {} - Number runs {}".format(flowcell.id, flowcell.number_runs)
    )
    logger.info(
        "Flowcell id: {} - Number barcodes {}".format(
            flowcell.id, flowcell.number_barcodes
        )
    )
    job_master.running = False
    job_master.save()
