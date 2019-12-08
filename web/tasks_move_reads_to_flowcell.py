from celery import task
from celery.utils.log import get_task_logger

from reads.models import (
    FastqRead,
    Flowcell,
    Run)

logger = get_task_logger(__name__)


@task()
def move_reads_to_flowcell(run_id, flowcell_id):
    """
    Moves the reads from flowcell A to flowcell B
    The split flowcell function calls this task
    """

    flowcell = Flowcell.objects.get(pk=flowcell_id)

    run = Run.objects.get(pk=run_id)

    reads = FastqRead.objects.filter(run=run)[:10000]

    reads_len = len(reads)

    while reads_len > 0:

        logger.info('Moving {} reads from run {} to flowcell {}'.format(
            reads_len,
            run.runid,
            flowcell.name
        ))

        print('Moving {} reads from run {} to flowcell {}'.format(
            reads_len,
            run.runid,
            flowcell.name
        ))

        for read in reads:
            read.flowcell = flowcell
            read.save()

        reads = FastqRead.objects.filter(run=run)[:10000]

        reads_len = len(reads)

    logger.info('Finished reads from run {} to flowcell {}'.format(
        run.runid,
        flowcell.name
    ))

    print('Finished reads from run {} to flowcell {}'.format(
        run.runid,
        flowcell.name
    ))
