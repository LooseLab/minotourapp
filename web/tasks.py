"""Contains code for some celery tasks
"""
from __future__ import absolute_import, unicode_literals

from celery.utils.log import get_task_logger

from alignment.tasks_alignment import run_minimap2_alignment
from artic.task_artic_alignment import (
    make_results_directory_artic,
    run_artic_command, run_artic_pipeline,
)
from metagenomics.new_centrifuge import run_centrifuge_pipeline
from metagenomics.sankey import calculate_sankey
from minknow_data.models import Flowcell
from minotourapp.celery import app
from minotourapp.utils import get_env_variable
from reads.models import (
    JobMaster,
)
from reads.tasks.redis_tasks_functions import harvest_reads
from reads.tasks.task_delete_flowcell import delete_flowcell
from reads.tasks.task_secure_artic_runs import secure_artic_runs
from reads.tasks.tasks_archive_flowcell import archive_flowcell
from reads.tasks.tasks_update_flowcell import update_flowcell_details
from readuntil.task_expected_benefit import calculate_expected_benefit_3dot0_final
from web.utils import fun

logger = get_task_logger(__name__)


@app.task
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

    # This fires a new task every 30 seconds to collect any uploaded reads and basically process them.
    harvest_reads.delay()
    if int(get_env_variable("MT_DESTROY_ARTIC_EVIDENCE")):
        secure_artic_runs.delay()

    # Create a list of all flowcells that have been active in the last 48 hours
    flowcell_list = [x for x in Flowcell.objects.all() if x.active()]
    # iterate them
    for flowcell in flowcell_list:
        if flowcell.archived:
            continue
        flowcell_job_list = JobMaster.objects.filter(flowcell=flowcell).filter(
            running=False, complete=False, paused=False
        )
        for flowcell_job in flowcell_job_list:
            if not flowcell_job.job_type.name == "Run Artic":
                flowcell_job.running = True
                flowcell_job.save()
            if flowcell_job.job_type.name == "Minimap2" and flowcell_job.from_database:
                run_minimap2_alignment.apply_async(args=(flowcell_job.id,), queue="minimap")
            if flowcell_job.job_type.name == "Metagenomics" and flowcell_job.from_database:
                run_centrifuge_pipeline.delay(flowcell_job.id)
            if flowcell_job.job_type.name == "UpdateFlowcellDetails":
                update_flowcell_details.delay(flowcell_job.id)
            if flowcell_job.job_type.name == "CalculateMetagenomicsSankey":
                run_sankey(flowcell_job.id)
            if flowcell_job.job_type.name == "Delete Flowcell":
                delete_flowcell.delay(flowcell_job.id)
            if flowcell_job.job_type.name == "Archive Flowcell":
                archive_flowcell.delay(flowcell_job.id)
            if flowcell_job.job_type.name == "ExpectedBenefit":
                calculate_expected_benefit_3dot0_final.delay(flowcell_job.id)
            if flowcell_job.job_type.name == "Track Artic Coverage" and flowcell_job.from_database:
                run_artic_pipeline.delay(flowcell_job.id)
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


@app.task(on_failure=fun)
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

