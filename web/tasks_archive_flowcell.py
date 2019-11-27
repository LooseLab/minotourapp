from celery import task
from celery.utils.log import get_task_logger
from jobs.models import JobMaster
from reads.models import FastqRead

logger = get_task_logger(__name__)

@task()
def archive_flowcell(task_id):
    """
    Asynchronous celery task to archive a flowcell, removing all fastq and rendering it permanently inactive. Data cannot be added to it.
    Parameters
    ----------
    task_id: int
        The primary key of the task (JobMaster)

    Returns
    -------
    None
    """

    task = JobMaster.objects.get(pk=task_id)

    flowcell = task.flowcell

    logger.info('Flowcell id: {} - Deleting flowcell {}'.format(flowcell.id, flowcell.name))

    # Get the first FastQRead object
    first_fastqread = FastqRead.objects.filter(flowcell=flowcell).order_by('id').first()
    # Get the last FastQRead object
    last_fastqread = FastqRead.objects.filter(flowcell=flowcell).order_by('id').last()
    # if we have both a first and last read PK
    if first_fastqread and last_fastqread:
        # Get the first fastqread primary key
        first_read = first_fastqread.id
        # Get the last fastqread primary key
        last_read = last_fastqread.id

        logger.info(
            'Flowcell id: {} - First and last fastqread id are {} {}'.format(flowcell.id, first_read, last_read))
        # Counter is the first read primary key
        counter = first_read

        # Whilst we still have reads left
        while counter <= last_read:
            counter = counter + 500

            logger.info('Flowcell id: {} - Deleting records with id < {}'.format(flowcell.id, counter))

            affected = FastqRead.objects.filter(flowcell=flowcell).filter(id__lt=counter).delete()

            logger.info('Flowcell id: {} - Deleted {} fastqread records'.format(flowcell.id, affected))

