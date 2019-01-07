from __future__ import absolute_import, unicode_literals

from celery import task
from celery.utils.log import get_task_logger
from django.db.models import Count, Max, Min, Sum, Avg

from reads.models import FastqRead, Flowcell, RunSummary

logger = get_task_logger(__name__)


@task()
def update_run_summary():

    for flowcell in Flowcell.objects.filter(is_active=True):

        for run in flowcell.runs.all():

            try:

                run_summary_list = RunSummary.objects.filter(run=run)

                if len(run_summary_list) > 0:

                    run_summary = run_summary_list[0]

                else:

                    run_summary = RunSummary.objects.create(

                        run=run,
                        read_count=None,
                        total_read_length=None,
                        max_read_length=None,
                        min_read_length=None,
                        avg_read_length=None,
                        last_read_start_time=None,
                        first_read_start_time=None,
                    )

            except RunSummary.DoesNotExist:

                run_summary = RunSummary.objects.create(

                    run=run,
                    read_count=None,
                    total_read_length=None,
                    max_read_length=None,
                    min_read_length=None,
                    avg_read_length=None,
                    last_read_start_time=None,
                    first_read_start_time=None,
                )

            if not run_summary.running:
                run_summary.running = True
                run_summary.save()

                reads = FastqRead.objects.filter(run=run).filter(id__gt=int(run_summary.last_read))[:50000]

                if reads.count() > 0:

                    result = reads.aggregate(
                        Count('sequence_length'),
                        Min('sequence_length'),
                        Max('sequence_length'),
                        Sum('sequence_length'),
                        #Avg('sequence_length'),
                        Min('start_time'),
                        Max('start_time'),
                        Max('id'),
                    )
                    if run_summary.read_count != None:
                        run_summary.read_count += result['sequence_length__count']
                        run_summary.total_read_length += result['sequence_length__sum']
                        run_summary.max_read_length = max(run_summary.max_read_length,result['sequence_length__max'])
                        run_summary.min_read_length = min(run_summary.min_read_length,result['sequence_length__min'])
                        #run_summary.avg_read_length = result['sequence_length__avg']
                        run_summary.first_read_start_time = min(run_summary.first_read_start_time,result['start_time__min'])
                        run_summary.last_read_start_time = max(run_summary.last_read_start_time,result['start_time__max'])
                    else:
                        run_summary.read_count = result['sequence_length__count']
                        run_summary.total_read_length = result['sequence_length__sum']
                        run_summary.max_read_length = result['sequence_length__max']
                        run_summary.min_read_length = result['sequence_length__min']
                        # run_summary.avg_read_length = result['sequence_length__avg']
                        run_summary.first_read_start_time =  result['start_time__min']
                        run_summary.last_read_start_time = result['start_time__max']
                    run_summary.last_read=result['id__max']
                    run_summary.avg_read_length=round(run_summary.total_read_length/run_summary.read_count)
                run_summary.running=False
                run_summary.save()

