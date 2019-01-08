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

            reads = FastqRead.objects.filter(run=run)

            if reads.count() > 0:

                result = reads.aggregate(

                    Count('sequence_length'),
                    Min('sequence_length'),
                    Max('sequence_length'),
                    Sum('sequence_length'),
                    Avg('sequence_length'),
                    Min('start_time'),
                    Max('start_time'),
                )

            else:

                result = {

                    'sequence_length__count': None,
                    'sequence_length__sum': None,
                    'sequence_length__max': None,
                    'sequence_length__min': None,
                    'sequence_length__avg': None,
                    'start_time__max': None,
                    'start_time__min': None,
                }

            try:

                run_summary_list = RunSummary.objects.filter(run=run)

                if len(run_summary_list) > 0:

                    run_summary = run_summary_list[0]

                else:

                    run_summary = RunSummary.objects.create(

                        run=run
                    )

            except RunSummary.DoesNotExist:

                run_summary = RunSummary.objects.create(

                    run=run
                )

            run_summary.read_count = result['sequence_length__count']
            run_summary.total_read_length = result['sequence_length__sum']
            run_summary.max_read_length = result['sequence_length__max']
            run_summary.min_read_length = result['sequence_length__min']
            run_summary.avg_read_length = result['sequence_length__avg']
            run_summary.first_read_start_time = result['start_time__min']
            run_summary.last_read_start_time = result['start_time__max']

            run_summary.save()

