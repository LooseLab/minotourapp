"""
Functionality to delete a previously run task with a celery task
"""
import datetime
import time
from pathlib import Path

from celery import task
from celery.utils.log import get_task_logger
from django.db import connection

from alignment.models import PafRoughCov, PafSummaryCov
from reads.models import JobMaster

# from readuntil.models import ExpectedBenefitChromosomes

logger = get_task_logger(__name__)


def get_or_create_results_directory(flowcell_id, task_id):
    """
    Make directory to store results of this task
    :param flowcell_id: The ID of the flowcell that produced these fastqs
    :param task_id: The ID of the Jobmaster that corresponds to this EB task
    :return: The path into this directory
    """
    results_dir = Path(f"readuntil/Temp_results/EB{flowcell_id}_{task_id}")
    if not results_dir.exists():
        Path.mkdir(results_dir)
    return results_dir


def delete_rows(first_row_object, rows_to_delete, job_master_object, reverse_lookup):
    """
    Delete as many rows from the database
    :param first_row_object: The django ORM object for the first row in the database to be deleted
    :param rows_to_delete: The number of rows to delete
    :param job_master_object: The JobMaster to delete these rows from, using a reverse lookup
    :param reverse_lookup: The key to lookup the objects keyed to the JobMaster
    :return: None
    """
    # Get the primary key of the first row
    first_row_object_pk = first_row_object.id
    # Let's delete the next thousand rows
    last_row_to_be_deleted = first_row_object_pk + rows_to_delete
    # delete them
    affected = getattr(job_master_object, reverse_lookup).filter(pk__lte=last_row_to_be_deleted).delete()

    logger.info('Flowcell id: {} - Deleted {} records'.format(job_master_object.flowcell.id, affected))


@task()
def delete_alignment_task(flowcell_job_id, restart=False):
    """
        The function called in monitor app to delete an alignment task,
        :param flowcell_job_id: The ID of the JobMaster database entry to be deleted
        :type flowcell_job_id: int
        :param restart: Restart the task after the deletion has occurred
        :type restart: bool
        :return:
    """
    # Get the flowcell jobMaster entry
    flowcell_job = JobMaster.objects.get(pk=flowcell_job_id)
    # set to complete so we don't generate new data
    flowcell_job.complete = True
    flowcell_job.save()
    # Get the flowcell
    flowcell = flowcell_job.flowcell
    delete = True
    delete_chunk_size = 100000
    paf_rough_cov_remaining = True
    while delete:
        if flowcell_job.job_type.name == "Minimap2" or flowcell_job.job_type.name == "Other":
            # If it's minimap2, we have tonnnes of alignment data, so delete it in batches of 1000
            # ######## Paf Rough cov ####### #
            if paf_rough_cov_remaining:
                # get the first entry in the database for this job master
                first_paf_rough_cov_row = PafRoughCov.objects.filter(job_master=flowcell_job)
                # if there isn't one we've deleted all the pafroughcov entries
                if not first_paf_rough_cov_row:
                    paf_rough_cov_remaining = False
                    logger.info(f"Deleted all Paf Rough Cov data. Job ID {flowcell_job.id}")
                else:
                    with connection.cursor() as cursor:
                        cursor.execute(
                            "DELETE FROM alignment_pafroughcov where flowcell_id=%s limit %s",
                            [flowcell.id, delete_chunk_size],
                        )
            if not paf_rough_cov_remaining:
                # Delete Paf Summary Cov
                PafSummaryCov.objects.filter(job_master=flowcell_job).delete()
                # if we are restarting the task, reset all the jobmaster fields to 0 and save them
                if restart:
                    flowcell_job.last_read = 0
                    flowcell_job.read_count = 0
                    flowcell_job.complete = False
                    flowcell_job.running = False
                    flowcell_job.save()
                    # We're restarting so set a new activity data so flowcell becomes active
                    flowcell.last_activity_date = datetime.datetime.now(datetime.timezone.utc)
                    flowcell.save()
                else:
                    finished = flowcell_job.delete()
                    logger.debug(f"Finished! {finished}")

                delete = False


@task()
def delete_metagenomics_task(flowcell_job_id, restart=False):
    """
    Delete the data from a metagenomics task, and then the database entry
    :param flowcell_job_id: The PK of the JobMaster database entry to be deleted
    :type flowcell_job_id: int
    :param restart: If this task is being restarted, we will clear all the data, delete the JobMaster from the
    database and reinsert it
    :type restart: bool
    :return:
    """
    # Get the flowcell jobMaster entry
    flowcell_job = JobMaster.objects.get(pk=flowcell_job_id)
    # Get the flowcell
    flowcell = flowcell_job.flowcell
    # Set complete to True so no new data is generated
    flowcell_job.complete = True
    flowcell_job.save()
    # If the task is running,it might be producing new data, so wait until it's finished
    running = flowcell_job.running
    seconds_counter = 0
    logger.info('Flowcell id: {} - Deleting alignment data from flowcell {}'.format(flowcell.id, flowcell.name))
    # keep looping
    delete = True
    centrifuge_output_remaining = True
    donut_data_remaining = True
    logger.info(f"Flowcell job type {flowcell_job.job_type.name}")
    if flowcell_job.job_type.name == "Metagenomics":

        # ## delete the alignment data from the target set validation ## #

        metagenomics_mapping_jobs = JobMaster.objects.filter(flowcell=flowcell, job_type__name="Other")
        # if there are mapping jobs
        if metagenomics_mapping_jobs:
            for mapping_job in metagenomics_mapping_jobs:
                logger.info(f"Deleting alignment data for {mapping_job.reference.name}")
                delete_alignment_task(mapping_job.id)
        else:
            logger.info("No mapping data found for this metagenomics run")
        # while we have data left
        while delete:
            if centrifuge_output_remaining:
                # Get the django ORM object representing the first row of the database entries to delete
                centrifuge_output_first_row = flowcell_job.centrifuge_output.first()
                # If there is no data left to delete
                if not centrifuge_output_first_row:
                    centrifuge_output_remaining = False
                else:
                    delete_rows(centrifuge_output_first_row, 1000, flowcell_job, "centrifuge_output")
            if donut_data_remaining:
                # Get the first row of the donut data database rows
                donut_data_first_row = flowcell_job.donut_data.first()
                # If there is no data left to delete
                if not donut_data_first_row:
                    donut_data_remaining = False
                else:

                    delete_rows(donut_data_first_row, 1000, flowcell_job, "donut_data")
            if not centrifuge_output_remaining and not donut_data_remaining:
                finished = flowcell_job.delete()
                logger.info(f"Finished deleting Task data for metagenomics {finished}")
                # if we are resetting the task reinsert a new fresh JobMaster
                if restart:
                    flowcell_job.last_read = 0
                    flowcell_job.read_count = 0
                    flowcell_job.iteration_count = 0
                    flowcell_job.complete = False
                    flowcell_job.save()
                    # We're restarting so set a new activity data so flowcell becomes active
                    flowcell.last_activity_date = datetime.datetime.now(datetime.timezone.utc)
                    flowcell.save()
                delete = False


@task()
def delete_expected_benefit_task(flowcell_job_id, restart=False):
    """
    Delete an expected benefit task along with the binary files it writes out
    :param flowcell_job_id: The Id of the JobMaster object to be deleted
    :param restart: Whether to begin the task again after deletion
    :return:
    """

    task = JobMaster.objects.get(pk=flowcell_job_id)

    if task is None:
        logger.error("Did not find JobMaster. Perhaps it has already been deleted?")
        return
    # Set to complete so no new data is produced
    running = task.running

    seconds_counter = 0

    while running:

        running = task.running
        time.sleep(1)
        seconds_counter += 1
        if seconds_counter == 140:
            task.running = False
            task.save()

    # task.save()
    # Get the flowcell
    flowcell = task.flowcell

    chromosomes = task.expectedBenefitTask.all()

    base_result_dir_path = get_or_create_results_directory(flowcell.id, task.id)

    for chromosome in chromosomes:
        chrom_key = chromosome.chromosome.line_name
        Path(f"{base_result_dir_path}/coverage_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
        Path(f"{base_result_dir_path}/counts_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
        Path(f"{base_result_dir_path}/mask_forward_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
        Path(f"{base_result_dir_path}/mask_reverse_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
        Path(f"{base_result_dir_path}/scores_forward_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
        Path(f"{base_result_dir_path}/scores_reverse_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
        Path(f"{base_result_dir_path}/cost_forward_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
        Path(f"{base_result_dir_path}/cost_reverse_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
        Path(f"{base_result_dir_path}/benefits_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
        Path(f"{base_result_dir_path}/fixed_benefits_forward{chrom_key}"
             f"_{flowcell.id}_{task.id}.dat").unlink()
        Path(f"{base_result_dir_path}/fixed_benefits_reverse{chrom_key}"
             f"_{flowcell.id}_{task.id}.dat").unlink()

    finished = task.delete()

    logger.info(f"Finished deleting Task data for expected benefit - {finished}")

    if restart:
        task.running = False
        task.complete = False
        task.last_read = 0
        task.read_count = 0
        task.iteration_count = 0
        
        task.save()
        flowcell.last_activity_date = datetime.datetime.now(datetime.timezone.utc)
        flowcell.save()
