from celery import task
from celery.utils.log import get_task_logger

from reads.models import (
    FastqRead,
    Flowcell,
    FlowcellChannelSummary,
    FlowcellStatisticBarcode,
    FlowcellHistogramSummary,
    FlowcellSummaryBarcode,
    Run)

import datetime

from reads.services import update_flowcell

logger = get_task_logger(__name__)


@task()
def move_reads_to_flowcell(run_id, flowcell_id, fromflowcell_id):
    """
    Moves the reads from flowcell A to flowcell B
    The split flowcell function calls this task
    """

    flowcell = Flowcell.objects.get(pk=flowcell_id)

    flowcellold = Flowcell.objects.get(pk=fromflowcell_id)

    run = Run.objects.get(pk=run_id)

    #reads = FastqRead.objects.filter(run=run)[:10000]
    chunksize = 100000

    totalreads = FastqRead.objects.filter(run=run).count()

    steps = int((totalreads - totalreads%chunksize)/chunksize + 1)

    for i in range(0, steps):

        modstart = i*chunksize

        modend = (i+1) * chunksize

        reads = FastqRead.objects.filter(run=run)[modstart:modend]

        reads_len = reads.count()

        logger.info('Moving {} reads from run {} to flowcell {}'.format(
            reads_len,
            run.runid,
            flowcell.name
        ))

        for read in reads:
            read.flowcell = flowcell
            read.save()

        update_flowcell.delay(reads)

    logger.info('Finished reads from run {} to flowcell {}'.format(
        run.runid,
        flowcell.name
    ))

    print('Finished reads from run {} to flowcell {}'.format(
        run.runid,
        flowcell.name
    ))

    ##Update tasks on old flowcell


    reset_flowcell_info(flowcellold)


def reset_flowcell_info(flowcell):
    # We need to delete the old objects for this flowcell
    FlowcellSummaryBarcode.objects.filter(flowcell=flowcell).delete()
    FlowcellStatisticBarcode.objects.filter(flowcell=flowcell).delete()
    FlowcellHistogramSummary.objects.filter(flowcell=flowcell).delete()
    FlowcellChannelSummary.objects.filter(flowcell=flowcell).delete()

    flowcell.start_time = datetime.datetime.now()
    flowcell.max_channel = 0
    flowcell.number_reads = 0
    flowcell.number_runs = 0
    flowcell.number_barcodes = 0
    flowcell.average_length = 0
    flowcell.total_read_length = 0
    flowcell.has_fastq = False
    flowcell.save()

    # recalculate old flowcell

    chunksize = 100000

    totalreads = FastqRead.objects.filter(flowcell=flowcell).count()

    steps = int((totalreads - totalreads % chunksize) / chunksize + 1)

    for i in range(0, steps):

        modstart = i * chunksize

        modend = (i + 1) * chunksize

        reads = FastqRead.objects.filter(flowcell=flowcell)[modstart:modend]

        reads_len = reads.count()

        logger.info('Moving {} reads from run {} to flowcell {}'.format(
            reads_len,
            run.runid,
            flowcell.name
        ))

        update_flowcell.delay(reads)

    logger.info('Finished reads from run {} to flowcell {}'.format(
        run.runid,
        flowcell.name
    ))



    #chancalcjob = flowcellold.flowcelljobs.get(job_type__name="ChanCalc")
    #chancalcjob.last_read=0
    #chancalcjob.read_count=0
    #chancalcjob.save()