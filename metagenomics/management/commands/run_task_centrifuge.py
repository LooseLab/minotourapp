"""
run_task_centrifuge.py
"""
import logging
import os

from django.core.management import BaseCommand, CommandError

from minknow_data.models import Flowcell
from reads.models import JobMaster
from web.tasks import run_centrifuge

log_folder = os.environ.get('MT_LOG_FOLDER')

logging.basicConfig(
    format='%(asctime)s %(module)s:%(levelname)s:%(thread)d:%(message)s',
    filename=os.path.join(log_folder, 'run_task_centrifuge.log'),
    level=os.environ.get('LOGLEVEL', 'INFO')
)

logger = logging.getLogger(__name__)

class Command(BaseCommand):

    help = 'Run task metagenomics'

    def handle(self, *args, **options):

        try:
            logger.info('Running metagenomics task')

            flowcell_list = Flowcell.objects.filter(is_active=True)

            for flowcell in flowcell_list:

                    flowcell_job_list = JobMaster.objects.filter(flowcell=flowcell).filter(running=False,
                                                                                           complete=False)

                    for flowcell_job in flowcell_job_list:

                        if flowcell_job.job_type.name == "Metagenomics":

                            logger.info("trying to run classification for flowcell {} {} {} ".format(
                                flowcell.id,
                                flowcell_job.id,
                                flowcell_job.last_read
                            ))

                            logger.info("starting metagenomics task")
                            run_centrifuge(flowcell_job.id)


        except Exception as e:
            logger.exception('Failed running task.')
            raise CommandError(repr(e))
