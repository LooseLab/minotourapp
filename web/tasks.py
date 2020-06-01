"""Contains code for some celery tasks
"""
from __future__ import absolute_import, unicode_literals

import redis
from celery import task
from celery.utils.log import get_task_logger

from alignment.tasks_alignment import run_minimap2_alignment
from artic.task_artic_alignment import (
    make_results_directory_artic,
    run_artic_command,
)
from centrifuge import centrifuge
from centrifuge.sankey import calculate_sankey
from minotourapp.utils import get_env_variable
from reads.models import (
    Run,
    Flowcell,
    JobMaster,
)
from reads.tasks.delete_flowcell_task import delete_flowcell
from reads.tasks.redis_tasks_functions import harvest_reads
from reads.tasks.update_flowcell_tasks import update_flowcell_details
from readuntil.task_expected_benefit import calculate_expected_benefit_3dot0_final

redis_instance = redis.StrictRedis(
    host="127.0.0.1", port=6379, db=0, decode_responses=True
)

logger = get_task_logger(__name__)


@task()
def run_monitor():
    """
    Run monitor is the task that is run every 30 seconds by Django celery beat -
    (minotourapp/settings.py for beat config)
    It iterates all the flowcells for flowcells that are active, getting all the attached JobMaster for jobs
    that are not running and not complete, and starts celery tasks for them.
    """

    logger.info("--------------------------------")
    logger.info("Running run_monitor celery task.")
    logger.info("--------------------------------")

    # This fires a new task every 90 seconds to collect any uploaded reads and basically process them.
    harvest_reads.delay()

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

            if flowcell_job.job_type.name == "Metagenomics":

                run_centrifuge.delay(flowcell_job.id)

            if flowcell_job.job_type.name == "UpdateFlowcellDetails":

                update_flowcell_details.delay(flowcell_job.id)

            if flowcell_job.job_type.name == "CalculateMetagenomicsSankey":

                run_sankey(flowcell_job.id)

            if flowcell_job.job_type.name == "Delete_Flowcell":

                delete_flowcell.delay(flowcell_job.id)

            if flowcell_job.job_type.name == "ExpectedBenefit":

                calculate_expected_benefit_3dot0_final.delay(flowcell_job.id)

            # if flowcell_job.job_type.name == "Track Artic Coverage":
            #
            #     run_artic_pipeline.delay(flowcell_job.id)

            if flowcell_job.job_type.name == "Run Artic":
                logger.info("Running the check for the Artic code.")
                # Get the maximum number of pipelines to have running
                max_artic_pipeline_instances = int(
                    get_env_variable("MT_ARTIC_MAX_PIPELINES")
                )

                logger.info(max_artic_pipeline_instances)
                # TODO replace with celery 4.0 magic queues etc
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

                    track_coverage_task_id = JobMaster.objects.get(
                        flowcell=flowcell_job.flowcell,
                        job_type__name="Track Artic Coverage",
                    ).id
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
def run_sankey(flowcell_job_id):
    """
    Calculate sankeys for a flowcell on demand of the user
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

