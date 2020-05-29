import logging
import os

from django.core.management.base import BaseCommand, CommandError

from alignment.tasks_alignment import run_minimap2_alignment_by_job_master
from reads.models import Flowcell, JobMaster

log_folder = os.environ.get('MT_LOG_FOLDER')

logging.basicConfig(
    format='%(asctime)s %(module)s:%(levelname)s:%(thread)d:%(message)s',
    filename=os.path.join(log_folder, 'run_task_minimap2.log'),
    level=os.environ.get('LOGLEVEL', 'INFO')
)

logger = logging.getLogger(__name__)


class Command(BaseCommand):

    help = 'Run minimap2 task'

    def handle(self, *args, **options):

        logger.info('Running minimap2 task')

        try:

            flowcell_list = Flowcell.objects.filter(is_active=True)

            for flowcell in flowcell_list:

                flowcell_job_list = JobMaster.objects.filter(flowcell=flowcell)

                for flowcell_job in flowcell_job_list:

                    if flowcell_job.job_type.name == "Minimap2":

                        run_minimap2_alignment_by_job_master(flowcell_job.id)

        except Exception as e:

            logger.exception('Failed running task.')
            raise CommandError(repr(e))
