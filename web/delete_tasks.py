"""
Functionality to delete a previously run task with a celery task
"""
from jobs.models import JobMaster
from celery import task
from celery.utils.log import get_task_logger
from pathlib import Path
import time
import datetime
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

    print('Flowcell id: {} - Deleted {} records'.format(job_master_object.flowcell.id, affected))


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
    
    # If the task is running,it might be producing new data, so wait until it's finished
    running = flowcell_job.running
    while running:
        running = flowcell_job.running
        time.sleep(3)
    
    # Get the flowcell
    flowcell = flowcell_job.flowcell

    logger.info('Flowcell id: {} - Deleting alignment data from flowcell {}'.format(flowcell.id, flowcell.name))
    # keep looping
    delete = True
    # Switch these to false when we have no task data left to delete
    pafstore_remaining = True

    paf_rough_cov_remaining = True

    while delete:

        if flowcell_job.job_type.name == "Minimap2":

            # If it's minimap2, we have tonnnes of alignment data, so delete it in batches of 1000

            # ######## Paf Rough cov ####### #
            if paf_rough_cov_remaining:
                # get the first entry in the database for this job master
                first_paf_rough_cov_row = flowcell_job.pafroughcov_list.first()
                # if there isn't one we've deleted all the pafroughcov entries
                if not first_paf_rough_cov_row:

                    paf_rough_cov_remaining = False

                else:
                    delete_rows(first_paf_rough_cov_row, 1000, flowcell_job, "pafroughcov_list")

            # ####### Paf Store ###### #
            if pafstore_remaining:
                # Get the first database entry
                first_pafstore_row = flowcell_job.pafstore_list.first()
                # If there is no data left
                if not first_pafstore_row:
                    print(f"finished deleting pafstore data")
                    pafstore_remaining = False

                else:

                    delete_rows(first_pafstore_row, 1000, flowcell_job, "paf_store_list")

            if not paf_rough_cov_remaining and not pafstore_remaining:
                # We've deleted all the super large data
                print(f"finished! {finshed}")

                finshed = flowcell_job.delete()

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
    while running:
        running = flowcell_job.running
        time.sleep(3)

    logger.info('Flowcell id: {} - Deleting alignment data from flowcell {}'.format(flowcell.id, flowcell.name))
    # keep looping
    delete = True

    centrifuge_output_remaining = True

    donut_data_remaining = True

    print(f"Flowcell job type {flowcell_job.job_type.name}")
    print(flowcell_job.job_type.name == "Metagenomics")

    if flowcell_job.job_type.name == "Metagenomics":

        # ## delete the alignment data from the target set validation ## #

        metagenomics_mapping_jobs = JobMaster.objects.filter(flowcell=flowcell, job_type__name="Other")
        print(metagenomics_mapping_jobs)
        # if there are mapping jobs
        if metagenomics_mapping_jobs:
            for mapping_job in metagenomics_mapping_jobs:
                print(f"deletin alignment data for {mapping_job.reference.name}")
                delete_alignment_task(mapping_job.id)

        else:
            print("no mapping data found for this metagenomics run")
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

            print()
            if not centrifuge_output_remaining and not donut_data_remaining:

                finished = flowcell_job.delete()

                print(f"Finished deleting Task data for metagenomics {finished}")

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
    # Set to complete so no new data is produced
    task.complete = True
    running = task.running
    # TODO do we wait until runnning is false so we don't have new data created or crash a running task?
    while running:
        running = task.running
        time.sleep(3)

    task.save()
    # Get the flowcell
    flowcell = task.flowcell

    # chromosomes = ExpectedBenefitChromosomes.objects.filter(task=task).values_list("chromosome__line_name", flat=True)

    base_result_dir_path = get_or_create_results_directory(flowcell.id, flowcell_job_id)

    # for chrom_key in chromosomes:
    #     # Write out the coverage
    #     Path(f"{base_result_dir_path}/coverage_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
    #     Path(f"{base_result_dir_path}/counts_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
    #     Path(f"{base_result_dir_path}/mask_forward_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
    #     Path(f"{base_result_dir_path}/mask_reverse_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
    #     Path(f"{base_result_dir_path}/scores_forward_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
    #     Path(f"{base_result_dir_path}/scores_reverse_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
    #     Path(f"{base_result_dir_path}/cost_forward_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
    #     Path(f"{base_result_dir_path}/cost_reverse_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
    #     Path(f"{base_result_dir_path}/benefits_{chrom_key}_{flowcell.id}_{task.id}.dat").unlink()
    #     Path(f"{base_result_dir_path}/fixed_benefits_forward{chrom_key}"
    #          f"_{flowcell.id}_{task.id}.dat").unlink()
    #     Path(f"{base_result_dir_path}/fixed_benefits_reverse{chrom_key}"
    #          f"_{flowcell.id}_{task.id}.dat").unlink()

    if restart:
        task.running = False
        task.complete = False
        task.last_read = 0
        task.read_count = 0
        task.iteration_count = 0
        
        task.save()
        flowcell.last_activity_date = datetime.datetime.now(datetime.timezone.utc)
        flowcell.save()
