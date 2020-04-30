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
    Flowcell,
    MinionRunInfo,
    JobMaster,
)
from web.tasks_chancalc import chancalc
from .tasks_alignment import run_minimap2_alignment
from pathlib import Path
from readuntil.task_expected_benefit import calculate_expected_benefit_3dot0_final


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
                    base_result_dir_path = make_results_directory_artic(
                        flowcell_job.flowcell.id, flowcell_job.id
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

    runs = Run.objects.all()

    for run in runs:

        if run.RunDetails.all().count():

            run.start_time = run.RunDetails.last().minKNOW_start_time
            origin = "Live data"
            run.save()
            print("Updating start_time for run {} from {}".format(run.runid, origin))


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