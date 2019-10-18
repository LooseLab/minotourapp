import logging
import os
import time

from django.core.management import BaseCommand, CommandError

from jobs.models import JobMaster
from reads.models import Flowcell
from web.tasks_chancalc import chancalc

log_folder = os.environ.get('MT_LOG_FOLDER')

logging.basicConfig(
    format='%(asctime)s %(module)s:%(levelname)s:%(thread)d:%(message)s',
    filename=os.path.join(log_folder, 'run_task_chancalc.log'),
    level=os.environ.get('LOGLEVEL', 'INFO')
)

logger = logging.getLogger(__name__)


class Command(BaseCommand):

    help = 'Run chancalc task'

    def handle(self, *args, **options):

        logger.info('Running chancalc task')

        try:
            flowcell_list = Flowcell.objects.filter(is_active=True)

            for flowcell in flowcell_list:

                flowcell_job_list = JobMaster.objects.filter(flowcell=flowcell)

                for flowcell_job in flowcell_job_list:

                    if flowcell_job.job_type.name == "ChanCalc":

                        chancalc(flowcell.id, flowcell_job.id, flowcell_job.last_read)

        except Exception as e:

            logger.exception('Failed running task.')
            raise CommandError(repr(e))
