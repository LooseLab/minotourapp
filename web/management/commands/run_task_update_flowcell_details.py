import logging
import os

from django.core.management import BaseCommand, CommandError

from reads.models import Flowcell, JobMaster
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

        logger.info('Running update flowcell details task')

        try:
            flowcell_list = Flowcell.objects.filter(is_active=True)

            for flowcell in flowcell_list:

                flowcell_job_list = JobMaster.objects.filter(flowcell=flowcell)

                for flowcell_job in flowcell_job_list:

                    if flowcell_job.job_type.name == "UpdateFlowcellDetails":

                        if not flowcell_job.running:

                            update_flowcell_details(flowcell.id, flowcell_job.id)



        except Exception as e:

            logger.exception('Failed running task.')
            raise CommandError(repr(e))
