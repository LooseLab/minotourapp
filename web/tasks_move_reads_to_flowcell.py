from celery import task
from celery.utils.log import get_task_logger

from reads.models import (
    FastqRead,
    Flowcell,
    Run)

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
    chunksize = 10000

    totalreads = FastqRead.objects.filter(run=run).count()

    steps = int((totalreads - totalreads%chunksize)/chunksize + 1)

    for i in range(0, steps):

        modstart = i*chunksize

        modend = (i+1) * chunksize

        reads = FastqRead.objects.filter(run=run)[modstart:modend]

        reads_len = FastqRead.objects.filter(run=run)[modstart:modend].count()

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


    logger.info('Finished reads from run {} to flowcell {}'.format(
        run.runid,
        flowcell.name
    ))

    print('Finished reads from run {} to flowcell {}'.format(
        run.runid,
        flowcell.name
    ))

    ##Update tasks on old flowcell

    chancalcjob = flowcellold.flowcelljobs.get(job_type__name="ChanCalc")

    chancalcjob.last_read=0
    chancalcjob.read_count=0
    chancalcjob.save()