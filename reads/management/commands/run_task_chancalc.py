from django.core.management import BaseCommand, CommandError

from devices.models import Flowcell
from jobs.models import JobMaster
from web.tasks import update_flowcell_details, processreads


class Command(BaseCommand):

    help = 'Run task chancalc'

    def handle(self, *args, **options):

        try:

            print('Running chancalc task')

            flowcell_list = Flowcell.objects.all()

            for flowcell in flowcell_list:

                flowcell_job_list = JobMaster.objects.filter(flowcell=flowcell)

                for flowcell_job in flowcell_job_list:

                    if flowcell_job.job_type.name == "ChanCalc":

                        processreads(flowcell.id, flowcell_job.id, flowcell_job.last_read)

        except Exception as e:

            raise CommandError(repr(e))