import datetime
import time

from celery.task import task
from celery.utils.log import get_task_logger
from django.db import connection

from minotourapp.utils import get_env_variable
from reads.models import JobMaster, FastqRead, Flowcell

logger = get_task_logger(__name__)

@task()
def archive_flowcell(flowcell_job_id):
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
    delete_chunk_size = 5000
    alignment_chunk_size = 7500
    logger.info(
        "Flowcell id: {} - Deleting flowcell {}".format(flowcell.id, flowcell.name)
    )
    JobMaster.objects.filter(flowcell=flowcell).update(complete=True)
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
            cursor.execute(
                "DELETE FROM reads_fastqread where flowcell_id=%s limit %s",
                [flowcell.id, delete_chunk_size],
            )
        logger.info(
            f"Flowcell id: {flowcell.id} - Deleted {delete_chunk_size} fastqread records in {time.time()-start_time}."
        )
        # TODO could this be rewritten into celery chunks?
        # time.sleep(3)
        archive_flowcell.apply_async(args=(flowcell_job_id,))

    else:
        #TODO sort out difference between archived flowcell so things to keep and just clearing a flowcell
        flowcell.archived = True
        flowcell.save()


@task()
def create_archive_tasks():
    """
    Create archive tasks for flowcells that are more than X days since last use, don't archive if set to -1
    Returns
    -------

    """
    time_until_inactive = get_env_variable("MT_TIME_UNTIL_ARCHIVE")
    if time_until_inactive == 0:
        return
    delta = datetime.timedelta(days=int(time_until_inactive))
    for flowcell in Flowcell.objects.filter(archived=False):
        if flowcell.last_activity_date < datetime.datetime.now(datetime.timezone.utc) - delta:
            JobMaster.objects.create(job_type_id=18, flowcell=flowcell)
            logger.info(f"Marking flowcell: {flowcell} for archiving.")
