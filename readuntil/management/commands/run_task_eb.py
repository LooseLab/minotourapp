import logging
import os

from django.core.management.base import BaseCommand, CommandError

from reads.models import Flowcell, JobType, JobMaster
from readuntil.task_expected_benefit import run_eb_by_job_master
from reference.models import ReferenceInfo

log_folder = os.environ.get('MT_LOG_FOLDER')

logging.basicConfig(
    format='%(asctime)s %(module)s:%(levelname)s:%(thread)d:%(message)s',
    filename=os.path.join(log_folder, 'run_task_minimap2.log'),
    level=os.environ.get('LOGLEVEL', 'INFO')
)

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = 'Run minimap2 task'

    def add_arguments(self, parser):
        """
        Add the arguments to the base command class
        :param parser:
        :return:
        """
        parser.add_argument(
            "job_master_id",
            help="The id of the Jobmaster for this EB task. If not extant, will be created.",
            nargs="+"
        )

    def handle(self, *args, **options):

        logger.info('Running minimap2 task')

        try:
            reference = ReferenceInfo.objects.get(name="Listeria_monocytogenes")
            flowcell = Flowcell.objects.all()[0]
            job_type = JobType.objects.get(name="ExpectedBenefit")
            job_master, created = JobMaster.objects.get_or_create(pk=options["job_master_id"],
                                                                  defaults={"flowcell": flowcell,
                                                                            "job_type": job_type,
                                                                            "reference": reference,
                                                                            "running": False

                                                                            })
            if created:
                print(f"Created Expected Benefit Jobmaster with PK {job_master.id}")

            else:
                print(f"Fetched Expected Benefit Jobmaster with PK {job_master.id}")

            print("Trying to run Expected Benefit task")
            run_eb_by_job_master(job_master.id)
            print("Scornfully begun task")

        except Exception as e:

            logger.exception('Failed running task.')
            raise CommandError(repr(e))
