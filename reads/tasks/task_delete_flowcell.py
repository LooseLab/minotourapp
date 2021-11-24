"""Delete a flowcell by first deleting it's reads, then it's database record
"""
import datetime
import time
from pathlib import Path
from shutil import rmtree

from celery.utils.log import get_task_logger
from django.db import connection

from alignment.models import PafSummaryCov
from minknow_data.models import Flowcell
from minotourapp.celery import app, MyTask
from minotourapp.utils import get_env_variable
from reads.models import (
    JobMaster,
    FastqRead,
    FastqFile,
    FlowcellChannelSummary,
    FlowcellHistogramSummary,
    FlowcellSummaryBarcode,
    FlowcellStatisticBarcode,
)
from web.delete_tasks import delete_alignment_task

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
        return 1
    else:
        # clear pngs from artic static dir
        rmtree(results_dir, onerror=on_delete_error)
        return 0


@app.task
def delete_flowcell(flowcell_id):
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
    flowcell_job = JobMaster.objects.filter(
        flowcell_id=flowcell_id, job_type_id=11
    ).last()
    # Get the flowcell
    flowcell = flowcell_job.flowcell
    delete_chunk_size = 5000
    logger.info(
        "Flowcell id: {} - Deleting flowcell {}".format(flowcell.id, flowcell.name)
    )
    JobMaster.objects.filter(flowcell=flowcell).update(complete=True)
    # Clear Artic data away
    for job in JobMaster.objects.filter(
        job_type__name="Track Artic Coverage", flowcell=flowcell
    ):
        clear_artic_data(job)
    # Get the last FastQRead object, wayyy faster than count
    last_fastq_read = FastqRead.objects.filter(flowcell=flowcell).last()
    # if we have both a first and last read PK
    if last_fastq_read:
        # Get the last fastqread primary key
        # Counter is the first 2read primary key
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
        delete_flowcell.apply_async(args=(flowcell_id,))
    else:
        if PafSummaryCov.objects.filter(job_master__in=flowcell.flowcell_jobs.all()).count():
            for fj in flowcell.flowcell_jobs.filter(job_type_id__in={4, 5}):
                delete_alignment_task(fj.id)
        affected = flowcell.delete()
        logger.info(f"Flowcell id: {flowcell.id} - Deleted {affected} flowcell.")


@app.task(base=MyTask)
def reset_flowcell(flowcell_pk, jm_pk):
    """
    Reset the basecalled data related aspects of a flowcell. Clears Fastq data, Fastq Files and all basecalled data summaries, and resets the flowcell values
    Parameters
    ----------
    flowcell_pk: int
        Primary key of the flowcell to be reset.
    jm_pk: int
        The primary key of the job master for this clear data
    Returns
    -------

    """
    flowcell = Flowcell.objects.get(pk=flowcell_pk)
    # delete fastq read records
    first_fastqread = FastqRead.objects.filter(flowcell=flowcell).first()
    delete_chunk_size = 5000
    if first_fastqread:
        with connection.cursor() as cursor:
            cursor.execute(
                "DELETE FROM reads_fastqread where flowcell_id=%s limit %s",
                [flowcell.id, delete_chunk_size],
            )
    else:
        flowcell.number_reads = 0
        flowcell.average_read_length = 0
        flowcell.total_read_length = 0
        flowcell.save()
        for run in flowcell.runs.all():
            if hasattr(run, "summary"):
                run.summary.delete()
        logger.info(f"Finished deleting fastq at {datetime.datetime.now()}")
        logger.info("Deleting Fastq file records...")
        affected = FastqFile.objects.filter(run__in=flowcell.runs.all()).delete()
        logger.info(f"Deleted: {affected}")
        logger.info("Deleting Basecalled data summaries...")
        affected = FlowcellChannelSummary.objects.filter(flowcell=flowcell).delete()
        logger.info(f"Deleted: {affected}")
        affected = FlowcellHistogramSummary.objects.filter(flowcell=flowcell).delete()
        logger.info(f"Deleted: {affected}")
        affected = FlowcellSummaryBarcode.objects.filter(flowcell=flowcell).delete()
        logger.info(f"Deleted: {affected}")
        affected = FlowcellStatisticBarcode.objects.filter(flowcell=flowcell).delete()
        logger.info(f"Deleted: {affected}")
        logger.info("Finished resetting flowcell!")
    left = FastqRead.objects.filter(flowcell=flowcell).count()
    logger.info(f"Time: {datetime.datetime.now()}, Reads Left: {left}")
    if left:
        reset_flowcell.apply_async(args=(flowcell_pk, jm_pk))
    else:
        JobMaster.objects.get(pk=jm_pk).delete()
