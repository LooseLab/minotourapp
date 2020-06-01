"""Delete a flowcell by first deleting it's reads, then it's database record
"""
import time

from celery.task import task
from celery.utils.log import get_task_logger
from django.db import connection

from alignment.models import PafStore
from reads.models import JobMaster, FastqRead

logger = get_task_logger(__name__)


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
