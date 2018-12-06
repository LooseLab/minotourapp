"""
run_task_sankey.py
"""

import logging
import os

from django.core.management import BaseCommand, CommandError

from jobs.models import JobMaster
from reads.models import Flowcell
from centrifuge.utils import calculate_sankey

log_folder = os.environ.get('MT_LOG_FOLDER')

logging.basicConfig(
    format='%(asctime)s %(module)s:%(levelname)s:%(thread)d:%(message)s',
    filename=os.path.join(log_folder, 'run_task_sankey.log'),
    level=os.environ.get('LOGLEVEL', 'INFO')
)

logger = logging.getLogger(__name__)

class Command(BaseCommand):

    help = 'Run task centrifuge'

    def handle(self, *args, **options):

        try:
            logger.info('Running sankey task')

            flowcell_list = Flowcell.objects.filter(is_active=True)

            for flowcell in flowcell_list:

                    flowcell_job_list = JobMaster.objects.filter(flowcell=flowcell).filter(running=False,
                                                                                           complete=False)

                    for flowcell_job in flowcell_job_list:

                        if flowcell_job.job_type.name == "CalculateSankey":

                            logger.info("trying to create sankey for flowcell {} {} {} ".format(
                                flowcell.id,
                                flowcell_job.id,
                                flowcell_job.last_read
                            ))

                            logger.info("starting Sankey task")
                            calculate_sankey(flowcell_job.id)


        except Exception as e:
            logger.exception('Failed running task.')
            raise CommandError(repr(e))
