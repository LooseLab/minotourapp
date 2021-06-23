import datetime

from celery.utils.log import get_task_logger

from artic.models import ArticBarcodeMetadata
from artic.task_artic_alignment import clear_unused_artic_files
from artic.utils import make_results_directory_artic
from minotourapp.celery import app
from minotourapp.utils import get_env_variable
from reads.models import JobMaster

logger = get_task_logger(__name__)


@app.task
def trigger_all_barcodes_after_run(artic_job_master):
    """
    Manually create run_artic_command job masters for the barcodes that never achieved sufficient coverage
    Parameters
    ----------
    artic_job_master: reads.models.JobMaster
        The artic job master on the inactive flowcell
    Returns
    -------

    """
    artic_barcodes = ArticBarcodeMetadata.objects.filter(
        job_master=artic_job_master, final_completion=False
    )
    for sub_cov_barcode in artic_barcodes:
        if sub_cov_barcode.barcode.name == "unclassified":
            continue
        sub_cov_barcode.final_completion = True
        sub_cov_barcode.save()
        logger.debug(f"Creating final artic run for {sub_cov_barcode} and then tidying up sensitive files...")
        # fixme - lazy here, we should really find a way to only create once from the sub cov barcode
        created, jm = JobMaster.objects.get_or_create(
            job_type_id=17,
            reference=artic_job_master.reference,
            barcode=sub_cov_barcode.barcode,
            flowcell=sub_cov_barcode.flowcell,
        )
        if not created:
            jm.complete = False
            jm.save()


@app.task
def secure_artic_runs():
    """
    When run monitor is celery beaten, go through flowcells with a run_artic_pipeline command and check they have had
     data uploaded to them in the last 6 hours.
    If they haven't, trigger any barcodes that haven't been run, to tidy up all the sensitive files.

    Returns
    -------

    """
    logger.info("Starting securing artic tasks for flowcells that haven't uploaded in 12 hours")
    jobs = JobMaster.objects.filter(job_type_id=16, complete=False)
    for artic_job in jobs:
        flowcell = artic_job.flowcell
        last_activity_date = flowcell.last_activity_date
        twelve_hours = datetime.timedelta(hours=int(get_env_variable("MT_ARTIC_TIME_UNTIL_CLEARING")))
        three_hours = datetime.timedelta(hours=3)
        active = (
                last_activity_date
                < datetime.datetime.now(datetime.timezone.utc) - twelve_hours
        )
        trigger_barcodes = (
            last_activity_date < datetime.datetime.now(datetime.timezone.utc) - three_hours
        )
        if trigger_barcodes:
            # so much unnecessary computing - need a way to mark it is final out side Artic Barcode Metadata
            trigger_all_barcodes_after_run(artic_job)

        if int(get_env_variable("MT_DESTROY_ARTIC_EVIDENCE")):
            # Not super happy with this, as it is affected by other things than read upload,
            # but can't think of an easy work around. If we aren't storing reads, we don't really update it apart from
            # when we upload a read batch
            # TOdo ideally we would add a run last activity time for this
            last_activity_date = flowcell.last_activity_date
            active = (
                last_activity_date
                > datetime.datetime.now(datetime.timezone.utc) - twelve_hours
            )
            if not active:
                results_dir = make_results_directory_artic(flowcell.id, artic_job.id)
                for barcode_name in ArticBarcodeMetadata.objects.filter(job_master=artic_job).values_list("barcode__name", flat=True):
                    clear_unused_artic_files(str(results_dir / barcode_name), barcode_name, flowcell.id)
                artic_job.complete = True
                artic_job.save()

    logger.info("Finished securing artic tasks")