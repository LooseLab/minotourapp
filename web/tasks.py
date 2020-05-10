from __future__ import absolute_import, unicode_literals
from celery import task
from celery.utils.log import get_task_logger
from django.db.models import Max, Count
from django.db import transaction
from artic.task_artic_alignment import make_results_directory_artic, run_artic_pipeline, run_artic_command
from centrifuge import centrifuge
from centrifuge.sankey import calculate_sankey
from minotourapp.utils import get_env_variable
from reads.models import (
    Barcode,
    FastqRead,
    Run,
    FlowcellSummaryBarcode,
    FlowcellStatisticBarcode,
    FlowcellHistogramSummary,
    Flowcell,
    MinionRunInfo,
    JobMaster,
    FlowcellChannelSummary,
)
from web.tasks_chancalc import chancalc
from .tasks_alignment import run_minimap2_alignment
from pathlib import Path
from readuntil.task_expected_benefit import calculate_expected_benefit_3dot0_final
from dateutil import parser

from reads.services import harvestreads

import redis
import numpy as np

redis_instance = redis.StrictRedis(host="127.0.0.1",
                                  port=6379, db=0,decode_responses=True)

logger = get_task_logger(__name__)


@task()
def run_monitor():
    """
    Run monitor is the task that is run every 30 seconds by Django celery beat -
    (minotourapp/settings.py for beat config)
    It iterates all the flowcells for flowcells that are active, getting all the attached JobMaster for jobs
    that are not running and not complete, and starts celery tasks for them.

    :return:
    """

    logger.info("--------------------------------")
    logger.info("Running run_monitor celery task.")
    logger.info("--------------------------------")

    # This fires a new task every 30 seconds to collect any uploaded reads and basically process them.
    harvestreads.delay()

    # Create a list of all flowcells that have been active in the last 48 hours
    flowcell_list = [x for x in Flowcell.objects.all() if x.active()]
    # iterate them
    for flowcell in flowcell_list:

        flowcell_job_list = JobMaster.objects.filter(flowcell=flowcell).filter(
            running=False, complete=False, paused=False
        )

        for flowcell_job in flowcell_job_list:

            if not flowcell_job.job_type.name == "Run Artic":
                flowcell_job.running = True
                flowcell_job.save()

            if flowcell_job.job_type.name == "Minimap2":

                run_minimap2_alignment.delay(flowcell_job.id,)

            if flowcell_job.job_type.name == "ChanCalc":

                chancalc.delay(flowcell.id, flowcell_job.id, flowcell_job.last_read)

            # if flowcell_job.job_type.name == "Assembly":
            #
            #     run_minimap2_assembly.delay(flowcell_job.id)

            if flowcell_job.job_type.name == "Metagenomics":

                run_centrifuge.delay(flowcell_job.id)

            if flowcell_job.job_type.name == "UpdateFlowcellDetails":

                update_flowcell_details.delay(flowcell_job.id)

            if flowcell_job.job_type.name == "CalculateMetagenomicsSankey":

                run_sankey(flowcell_job.id)

            if flowcell_job.job_type.name == "Delete_Flowcell":

                run_delete_flowcell.delay(flowcell_job.id)

            if flowcell_job.job_type.name == "ExpectedBenefit":

                calculate_expected_benefit_3dot0_final.delay(flowcell_job.id)

            if flowcell_job.job_type.name == "Track Artic Coverage":

                run_artic_pipeline.delay(flowcell_job.id)

            if flowcell_job.job_type.name == "Run Artic":
                logger.info("Running the check for the Artic code.")
                # Get the maximum number of pipelines to have running
                max_artic_pipeline_instances = int(
                    get_env_variable("MT_ARTIC_MAX_PIPELINES")
                )

                logger.info(max_artic_pipeline_instances)
                # TODO replace with celery 4.0 magic
                running_artic_jobs = JobMaster.objects.filter(
                    job_type=flowcell_job.job_type, running=True
                ).count()

                logger.info(
                    f" Running this many artic jobs currently {running_artic_jobs}"
                )

                if running_artic_jobs < max_artic_pipeline_instances:

                    flowcell_job.running = True

                    flowcell_job.save()

                    logger.info("Running the task for the Artic code.")

                    barcode_name = flowcell_job.barcode.name
                    # generate the base results path

                    track_coverage_task_id = JobMaster.objects.get(flowcell=flowcell_job.flowcell
                                                                   , job_type__name="Track Artic Coverage").id
                    base_result_dir_path = make_results_directory_artic(
                        flowcell_job.flowcell.id, track_coverage_task_id
                    )

                    # run the command task.
                    run_artic_command.delay(
                        str(base_result_dir_path), barcode_name, flowcell_job.id
                    )

                    logger.info("Added Artic to the redis cli queue.")

                elif running_artic_jobs >= max_artic_pipeline_instances:
                    logger.warning(
                        f"Sadly we do not have enough available compute power to run more than {max_artic_pipeline_instances} Artic pipelines simultaneously"
                    )

                else:
                    logger.error("¯\_(ツ)_/¯")

@task()
def run_delete_flowcell(flowcell_job_id):
    """
    The function called in monitorapp to delete a flowcell
    :param flowcell_job_id: The primary key of the JobMaster database entry
    :return:
    """

    # Get the flowcell jobMaster entry
    flowcell_job = JobMaster.objects.get(pk=flowcell_job_id)
    # Get the flowcell
    flowcell = flowcell_job.flowcell

    logger.info(
        "Flowcell id: {} - Deleting flowcell {}".format(flowcell.id, flowcell.name)
    )

    # Get the first FastQRead object
    first_fastqread = FastqRead.objects.filter(flowcell=flowcell).order_by("id").first()
    # Get the last FastQRead object
    last_fastqread = FastqRead.objects.filter(flowcell=flowcell).order_by("id").last()
    # if we have both a first and last read PK
    if first_fastqread and last_fastqread:
        # Get the first fastqread primary key
        first_read = first_fastqread.id
        # Get the last fastqread primary key
        last_read = last_fastqread.id

        logger.info(
            "Flowcell id: {} - First and last fastqread id are {} {}".format(
                flowcell.id, first_read, last_read
            )
        )
        # Counter is the first read primary key
        counter = first_read

        # Whilst we still have reads left
        while counter <= last_read:
            counter = counter + 5000

            logger.info(
                "Flowcell id: {} - Deleting records with id < {}".format(
                    flowcell.id, counter
                )
            )

            affected = (
                FastqRead.objects.filter(flowcell=flowcell)
                .filter(id__lt=counter)
                .delete()
            )

            logger.info(
                "Flowcell id: {} - Deleted {} fastqread records".format(
                    flowcell.id, affected
                )
            )

    affected = Flowcell.objects.get(pk=flowcell.id).delete()

    logger.info("Flowcell id: {} - Deleted {} records".format(flowcell.id, affected))


@task()
def run_sankey(flowcell_job_id):
    """
    Calculate sankeys for a flowcell upon demand of the user
    :param flowcell_job_id: The pk id of this task
    :author: Rory
    :return:
    """
    job_master = JobMaster.objects.get(pk=flowcell_job_id)
    logger.info(
        "Flowcell id: {} - Starting Sankey calculation task".format(
            job_master.flowcell.id
        )
    )
    calculate_sankey(flowcell_job_id)


@task()
def run_centrifuge(flowcell_job_id):

    job_master = JobMaster.objects.get(pk=flowcell_job_id)

    logger.info(
        "Flowcell id: {} - Starting centrifuge task".format(job_master.flowcell.id)
    )

    centrifuge.run_centrifuge(flowcell_job_id)


@task()
def delete_runs():
    """
    Delete runs that have been marked for deletion. Called by Celery beat.
    Returns
    -------

    """
    Run.objects.filter(to_delete=True).delete()


@task()
def update_run_start_time():
    """
    This method update the field start_time of run based on the live data or
    the header of the first fastq read
    """
    # TODO could this be improved?

    runs = Run.objects.all()

    for run in runs:

        if run.RunDetails.all().count():

            run.start_time = run.RunDetails.last().minKNOW_start_time
            origin = "Live data"
            run.save()
            logger.info(f"Updating start_time for run {run.runid} from {origin}")


def scan_keys(r, pattern):
    "Returns a list of all the keys matching a given pattern in a redis index"
    result = []
    cur, keys = r.scan(cursor=0, match=pattern, count=2)
    result.extend(keys)
    while cur != 0:
        cur, keys = r.scan(cursor=cur, match=pattern, count=2)
        result.extend(keys)
    return result


def gd_key(r, key):
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


def hgd_key(r, key):
    """
    Safely delete a redis hash (Dictionary)
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
    p.hgetall(key)
    p.delete(key)
    returnvalue = p.execute()[0]
    return returnvalue

@task(serializer="pickle")
def fsumb(flowcell):
    """
    Flowcell summary barcode, get all keys related to the flowcell summary barcode values of this flowcell
    from redis and update them
    Parameters
    ----------
    flowcell

    Returns
    -------

    """
    keys = redis_instance.scan_iter("{}_flowcellSummaryBarcode_*".format(flowcell.id))
    for key in keys:
        result = hgd_key(redis_instance, key)
        barcode_name = result['barcode_name']
        rejection_status = result['rejection_status']
        type_name = result['type_name']
        status = result['status']
        sequence_length_sum = int(result['total_length'])
        quality_average_sum = int(result['quality_sum'])
        read_count = int(result['read_count'])
        sequence_length_max = int(result["max_length"])
        sequence_length_min = int(result["min_length"])
        new_channel_list = list(result["channel_presence"])

        flowcellSummaryBarcode, created = FlowcellSummaryBarcode.objects.get_or_create(
            flowcell=flowcell,
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
        if flowcellSummaryBarcode.min_length == 0 or flowcellSummaryBarcode.min_length > sequence_length_min:
            flowcellSummaryBarcode.min_length = sequence_length_min

        channel_list = list(flowcellSummaryBarcode.channel_presence)
        fusedlist = np.bitwise_or(np.array(channel_list, dtype=int), np.array(new_channel_list, dtype=int))

        flowcellSummaryBarcode.channel_presence = ("").join(list(fusedlist.astype(str)))
        flowcellSummaryBarcode.channel_count = np.count_nonzero(fusedlist == 1)

        flowcellSummaryBarcode.save()


@task(serializer="pickle")
def fhs(flowcell):
    """
    Flowcell histogram summary. Get all keys related to the flowcell histogram values of this flowcell
    from redis and update them
    Parameters
    ----------
    flowcell

    Returns
    -------

    """
    keys = redis_instance.scan_iter("{}_flowcellHistogramSummary_*".format(flowcell.id))
    for key in keys:
        result = hgd_key(redis_instance, key)
        bin_index = int(result['bin_index'])
        barcode_name = result['barcode_name']
        rejection_status = result['rejection_status']
        read_type_name = result['read_type_name']
        status = result['status']
        sequence_length_sum = int(result['read_length'])
        read_count = int(result['read_count'])

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

@task(serializer="pickle")
def fsb(flowcell):
    """
    Flowcell statistic Barcode. Get all keys related to the flowcell summary barcode values of this flowcell
    from redis and update them
    Parameters
    ----------
    flowcell

    Returns
    -------

    """
    keys = redis_instance.scan_iter("{}_flowcellStatisticBarcode_*".format(flowcell.id))
    for key in keys:
        result = hgd_key(redis_instance, key)
        start_time = parser.parse(result['sample_time'])
        barcode_name = result['barcode_name']
        rejection_status = result['rejection_status']
        type_name = result['type_name']
        status = result['status']
        sequence_length_sum = int(result['total_length'])
        quality_average_sum = int(result['quality_sum'])
        read_count = int(result['read_count'])
        sequence_length_max = int(result["max_length"])
        sequence_length_min = int(result["min_length"])
        new_channel_list = list(result["channel_presence"])

        flowcellStatisticBarcode, created = FlowcellStatisticBarcode.objects.get_or_create(
            flowcell=flowcell,
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
        fusedlist = np.bitwise_or(np.array(channel_list, dtype=int), np.array(new_channel_list, dtype=int))

        flowcellStatisticBarcode.channel_presence = ("").join(list(fusedlist.astype(str)))
        flowcellStatisticBarcode.channel_count = np.count_nonzero(fusedlist == 1)

        flowcellStatisticBarcode.save()

@task()
def update_flowcell_details(job_master_id):
    """
    This task updates the flowcell details (number of runs, number of reads, sample name)
    using the MinIONRunStatus records if they are available, otherwise reading from the
    fastq file summaries
    """
    # TODO
    """
    Refactoring this task such that:
    1. it only runs when new reads have appeared.
    2. it only processes reads since the last one was processed.
    3. Need to look at the reset on this task - it doesn't seem to work.                    

    """

    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = True
    job_master.save()

    with transaction.atomic():
        flowcell = job_master.flowcell

        logger.info('Flowcell id: {} - Updating details of flowcell {}'.format(flowcell.id, flowcell.name))

        #
        # Get the first MinIONRunStatus for a particular flowcell - but we want to ignore platform QCs
        #
        minion_run_status_first = MinionRunInfo.objects.filter(run_id__flowcell=flowcell).exclude(
            experiment_type="platform_qc").order_by('minKNOW_start_time') \
            .first()

        #
        # If the MinIONRunStatus exists, than update start time
        #
        if minion_run_status_first:
            logger.info('Flowcell id: {} - There is at least one MinIONRunStatus'.format(flowcell.id))

            flowcell.start_time = minion_run_status_first.minKNOW_start_time

            logger.info('Flowcell id: {} - Setting start_time to {}'.format(flowcell.id, flowcell.start_time))
        #
        # Alternatively check the reads to find the earliest read?
        #
        #else:
        #    first_read_entry =

        #
        # Define flowcell's sample_name
        #
        minion_run_status_list = MinionRunInfo.objects.filter(run_id__flowcell=flowcell).exclude(
            experiment_type="platform_qc")

        if minion_run_status_list.count() > 0:

            for minion_run_status in minion_run_status_list:

                if minion_run_status.minKNOW_sample_name and minion_run_status.minKNOW_sample_name != 'undefined':
                    flowcell.sample_name = minion_run_status_first.minKNOW_sample_name

                    logger.info(
                        'Flowcell id: {} - Setting sample_name to {} - data from MinIONRunStatus.minKNOW_sample_name'.format(
                            flowcell.id, flowcell.sample_name))

                    break

        else:

            for run in Run.objects.filter(flowcell=flowcell).exclude(name="mux scan"):

                if run.name and run.name != 'undefined':
                    flowcell.sample_name = run.name

                    logger.info('Flowcell id: {} - Setting sample_name to {} - data from Run.name'.format(flowcell.id,
                                                                                                          flowcell.sample_name))

                    break

        #
        # Get number of fastqreads
        #
        number_reads = job_master.read_count

        last_read = job_master.last_read

        # Todo: If the flowcell has no reads, it may have been reset and so we may need to recount it based on previous issue.

        #if redis_instance.get("{}_has_fastq".format(flowcell.id)):
        #    flowcell.has_fastq = True

        #flowcell.size = redis_instance.get("{}_size".format(flowcell.id))

        #redis_instance.delete("{}_size".format(flowcell.id))

        #max_channel = redis_instance.get("{}_max_channel".format(flowcell.id))
        #if max_channel:
        #    if max_channel > flowcell.max_channel:
        #        flowcell.max_channel = max_channel
        #    redis_instance.delete("{}_max_channel".format(flowcell.id))

        ### Get flowcell max_channel number.
        has_fastq = gd_key(redis_instance, "{}_has_fastq".format(flowcell.id))
        if has_fastq:
            flowcell.has_fastq = has_fastq

            total_read_length = gd_key(redis_instance,"{}_total_read_length".format(flowcell.id))
            if total_read_length:
                flowcell.total_read_length += int(total_read_length)

            number_reads = gd_key(redis_instance, "{}_number_reads".format(flowcell.id))
            if number_reads:
                flowcell.number_reads += int(number_reads)

            if number_reads and total_read_length:
                flowcell.average_read_length = int(total_read_length)/int(number_reads)

            max_channel = gd_key(redis_instance,"{}_max_channel".format(flowcell.id))
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


            ### Now to try and update all the channel values that we need to update...

            for channel in range(1,flowcell.max_channel+1):
                read_count = gd_key(redis_instance,"{}_{}_read_count".format(flowcell.id,channel))
                read_length = gd_key(redis_instance, "{}_{}_read_length".format(flowcell.id, channel))
                if read_count or read_length:
                    flowcellChannelSummary, created = FlowcellChannelSummary.objects.get_or_create(
                        flowcell=flowcell, channel=channel
                    )
                    if read_count:
                        flowcellChannelSummary.read_count += int(read_count)
                    if read_length:
                        flowcellChannelSummary.read_length += int(read_length)
                    flowcellChannelSummary.save()


            ### Now try and update all the flowcellsummarystatistics

            fsb.delay(flowcell)

            fhs.delay(flowcell)

            fsumb.delay(flowcell)
            



        ### This query is slow - and it should be fast.

        #
        # # Update flowcell size
        # #
        ##max_channel = FastqRead.objects.filter(
        ##    flowcell=flowcell,
        ##    id__gt=int(job_master.last_read)
        ##).aggregate(result=Max('channel'), last_read=Max('id'),read_number=Count('read_id'))

        ##if max_channel["read_number"] is not None:

        ##    number_reads += max_channel["read_number"]

        ##    flowcell.has_fastq = True

        ##job_master.read_count = number_reads

        ##if max_channel['last_read'] is not None:

        ##    job_master.last_read = max_channel['last_read']

        # job_master.save()
        #
        # Get the job_master chancalc for this flowcell
        #

        ##job_master_list = JobMaster.objects.filter(flowcell=flowcell, job_type__name='Chancalc')

        ##number_reads_processed = flowcell.number_reads_processed

        ##if job_master_list.count() > 0:

        ##    number_reads_processed = job_master_list[0].read_count

        #
        # Get the FlowcellSummaryBarcodes for a particular flowcell and for barcode_name "All reads"
        #
        ##flowcell_summary_list = FlowcellSummaryBarcode.objects.filter(flowcell=flowcell).filter(barcode_name='All reads')

        ##average_read_length = 0

        ##total_read_length = 0

        ##logger.info('Flowcell id: {} - There is/are {} FlowcellSummaryBarcode records'.format(flowcell.id, len(flowcell_summary_list)))

        ##for flowcell_summary in flowcell_summary_list:

        ##    total_read_length += flowcell_summary.total_length

        ##if number_reads > 0:

        ##    average_read_length = total_read_length / number_reads

        ##logger.info('Flowcell id: {} - Total read length {}'.format(flowcell.id, total_read_length))
        ##logger.info('Flowcell id: {} - Number reads {}'.format(flowcell.id, number_reads))
        ##logger.info('Flowcell id: {} - Number reads processed {}'.format(flowcell.id, number_reads_processed))
        ##logger.info('Flowcell id: {} - Average read length {}'.format(flowcell.id, average_read_length))

        ##flowcell.average_read_length = average_read_length
        ##flowcell.total_read_length = total_read_length
        ##flowcell.number_reads = number_reads
        ##flowcell.number_reads_processed = number_reads_processed

        flowcell.number_runs = Run.objects.filter(flowcell=flowcell).count()

        barcode_count = len(Barcode.objects.filter(run_id__flowcell=flowcell).values('name').distinct()) - 1

        flowcell.number_barcodes = barcode_count

        flowcell.save()

        logger.info('Flowcell id: {} - Number runs {}'.format(flowcell.id, flowcell.number_runs))
        logger.info('Flowcell id: {} - Number barcodes {}'.format(flowcell.id, flowcell.number_barcodes))

        # #
        # # Update flowcell size
        # #
        # max_channel = FastqRead.objects.filter(
        #     flowcell=flowcell,
        #     id__gt=int(job_master.last_read)
        # ).aggregate(result=Max('channel'), last_read=Max('id'),read_number=Count('read_id'))

        ##if max_channel['result'] is not None and max_channel['result'] > flowcell.max_channel:

        ##    flowcell.max_channel = max_channel['result']

        ##    if max_channel['result']:

        ##        if max_channel['result'] > 512:

        ##            flowcell.size = 3000

        ##        elif max_channel['result'] > 126:

        ##            flowcell.size = 512

        ##        else:

        ##            flowcell.size = 126

        ##    else:

        ##        flowcell.size = 512

        ##flowcell.save()

        # job_master = JobMaster.objects.get(pk=job_master_id)

        job_master.running = False

        job_master.save()