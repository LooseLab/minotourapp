"""
Functionality to delete a previously run task with a celery task
"""
from jobs.models import JobMaster
from celery import task
from celery.utils.log import get_task_logger


logger = get_task_logger(__name__)


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

    print('Flowcell id: {} - Deleted {} fastqread records'.format(job_master_object.flowcell.id, affected))

@task()
def delete_alignment_task(flowcell_job_id):
    """
        The function called in monitor app to delete an alignment task,
        :param flowcell_job_id: The ID of the JobMaster database entry to be deleted
        :return:
    """

    # Get the flowcell jobMaster entry
    flowcell_job = JobMaster.objects.get(pk=flowcell_job_id)
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
                finshed = flowcell_job.delete()

                print(f"finiehsed! {finshed}")

                delete = False


@task()
def delete_metagenomics_task(flowcell_job_id):
    """
    Delete the data from a metagenomics task, and then the database entry
    :param flowcell_job_id: The PK of the JobMaster database entry to be deleted
    :return:
    """

    # Get the flowcell jobMaster entry
    flowcell_job = JobMaster.objects.get(pk=flowcell_job_id)
    # Get the flowcell
    flowcell = flowcell_job.flowcell

    logger.info('Flowcell id: {} - Deleting alignment data from flowcell {}'.format(flowcell.id, flowcell.name))
    # keep looping
    delete = True

    centrifuge_output_remaining = True

    donut_data_remaining = True

    if flowcell_job.job_type.name == "Metagenomics":

        # ## delete the alignment data from the target set validation ## #

        delete_alignment_task(flowcell_job_id)
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

                print(f"Finished deleting Task data for metagenomics {finished}")

                delete = False

def delete_expected_benefit_task(flowcell_job_id):
    """
    Delete an expected benefit task along with the binary files it writes out
    :param flowcell_job_id: The Id of the JobMaster object to be deleted
    :return:
    """