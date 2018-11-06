import logging
import os

from django.core.management import BaseCommand, CommandError

from web.tasks import update_flowcell_details


log_folder = os.environ.get('MT_LOG_FOLDER')

logging.basicConfig(
    format='%(asctime)s %(module)s:%(levelname)s:%(thread)d:%(message)s',
    filename=os.path.join(log_folder, 'run_task_update_flowcell_details.log'),
    level=os.environ.get('LOGLEVEL', 'INFO')
)

logger = logging.getLogger(__name__)


class Command(BaseCommand):

    help = 'Run update flowcell details task'

    def handle(self, *args, **options):

        try:

            logger.info('Running update flowcell details task')

            update_flowcell_details()

        except Exception as e:

            logger.exception('Failed running task.')
            raise CommandError(repr(e))
