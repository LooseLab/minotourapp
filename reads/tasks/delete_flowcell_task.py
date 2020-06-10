"""Delete a flowcell by first deleting it's reads, then it's database record
"""
import time
from pathlib import Path
from shutil import rmtree

from celery.task import task
from celery.utils.log import get_task_logger
from django.db import connection

from alignment.models import PafStore
from minotourapp.settings import BASE_DIR
from minotourapp.utils import get_env_variable
from reads.models import JobMaster, FastqRead

logger = get_task_logger(__name__)

def on_delete_error(function, path, excinfo):
    """
    Function provided to shutils rmtree on error catch
    Parameters
    ----------
    function
        The function that raised the exception
    path
        Path name passed to function
    excinfo
        The exception information returned by sys.execinfo

    Returns
    -------

    """
    logger.error(excinfo)
    return 1


def clear_artic_data(job_master):
    """
    Clear the artic files from the system drive
    Parameters
    ----------
    job_master: reads.models.JobMaster
        The job master ORM object of the track artic job
    Returns
    -------
    exit_code: int
        0 if successful, 1 if not
    """
    environmental_results_directory = get_env_variable("MT_ARTIC_RESULTS_DIR")
    results_dir = Path(
        f"{environmental_results_directory}/artic/Temp_results/{job_master.flowcell.id}_{job_master.id}_artic"
    )

    if not results_dir.exists():
        exit_code = 1
    else:
        # clear pngs from artic static dir
        rmtree(Path(f"{BASE_DIR}/artic/static/artic/{job_master.flowcell.id}_{job_master.id}_artic/"), onerror=on_delete_error)
        rmtree(results_dir, onerror=on_delete_error)
        return 0


@task()
def delete_flowcell(flowcell_job_id):
    """
    Delete all reads from a flowcell, then delete the flowcell itself.
    Parameters
    ----------
    flowcell_job_id: int
        The PK of the job_master database record of this task

    Returns
    -------

    """
    # Get the flowcell jobMaster entry
    start_time = time.time()
    flowcell_job = JobMaster.objects.get(pk=flowcell_job_id)
    # Get the flowcell
    flowcell = flowcell_job.flowcell
    delete_chunk_size = 100000
    alignment_chunk_size = 150000
    logger.info(
        "Flowcell id: {} - Deleting flowcell {}".format(flowcell.id, flowcell.name)
    )

    JobMaster.objects.filter(flowcell=flowcell).update(complete=True)
    # Clear Artic data away
    for job in JobMaster.objects.filter(job_type__name="Track Artic Coverage", flowcell=flowcell):
        clear_artic_data(job)

    # Get the last FastQRead object, wayyy faster than count
    last_fastq_read = FastqRead.objects.filter(flowcell=flowcell).last()
    # if we have both a first and last read PK
    if last_fastq_read:
        # Get the last fastqread primary key
        last_read = last_fastq_read.id
        # Counter is the first read primary key
        # Whilst we still have reads left

        # counter = counter + 20000
        logger.info(
            f"Flowcell id: {flowcell.id} - Deleting {delete_chunk_size} records."
        )
        with connection.cursor() as cursor:
            if PafStore.objects.filter(job_master__flowcell=flowcell).first():
                # TODO this is slow - we want to store flowcell ID on pafstore, there 's already a field, just not used
                cursor.execute(
                    "DELETE alignment_pafstore FROM alignment_pafstore INNER JOIN reads_jobmaster ON alignment_pafstore.job_master_id = reads_jobmaster.id WHERE reads_jobmaster.flowcell_id=%s and alignment_pafstore.id < %s",
                    [flowcell.id, PafStore.objects.filter(job_master__flowcell=flowcell).first().id + alignment_chunk_size],
                )
                logger.info(f"Deleted alignment entries.")
            cursor.execute(
                "DELETE FROM reads_fastqread where flowcell_id=%s limit %s",
                [flowcell.id, delete_chunk_size],
            )
        logger.info(
            f"Flowcell id: {flowcell.id} - Deleted {delete_chunk_size} fastqread records in {time.time()-start_time}."
        )
        # TODO could this be rewritten into celery chunks?
        time.sleep(3)
        delete_flowcell.apply_async(args=(flowcell_job_id,))
    else:
        affected = flowcell.delete()
        logger.info(f"Flowcell id: {flowcell.id} - Deleted {affected} flowcell.")
