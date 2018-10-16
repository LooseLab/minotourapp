from django.core.management.base import BaseCommand, CommandError

from jobs.models import JobMaster
from web.tasks_alignment import run_minimap2_alignment_by_job_master


class Command(BaseCommand):

    help = 'Run minimap2 task for a specific flowcell'

    def add_arguments(self, parser):
        parser.add_argument('task_id', type=int)

    def handle(self, *args, **options):

        print('Recovering information about the task')

        task = JobMaster.objects.get(pk=int(options['task_id']))

        print('---> task id: {}'.format(task.id))
        print('---> flowcell id: {}'.format(task.flowcell.id))
        print('---> reference id: {}'.format(task.reference.id))
        print('---> last read: {}'.format(task.last_read))

        run_minimap2_alignment_by_job_master(task.id)
