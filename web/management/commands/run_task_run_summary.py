import logging
import os

from django.core.management import BaseCommand, CommandError

from reads.tasks.tasks_update_run_summary import update_run_summary

log_folder = os.environ.get('MT_LOG_FOLDER')

logging.basicConfig(
    format='%(asctime)s %(module)s:%(levelname)s:%(thread)d:%(message)s',
    filename=os.path.join(log_folder, 'run_task_run_summary.log'),
    level=os.environ.get('LOGLEVEL', 'INFO')
)

logger = logging.getLogger(__name__)


class Command(BaseCommand):

    help = 'Run run_summary task'

    def handle(self, *args, **options):

        logger.info('Running run_summary task')

        try:

            update_run_summary()

        except Exception as e:

            logger.exception('Failed running task.')
            raise CommandError(repr(e))
