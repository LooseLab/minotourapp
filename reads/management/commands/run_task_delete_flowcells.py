from django.core.management import BaseCommand, CommandError

from reads.models import Flowcell
from jobs.models import JobMaster
from web.tasks import update_flowcell_details, processreads


class Command(BaseCommand):

    help = 'Run task delete flowcells'

    def handle(self, *args, **options):

        print('Run task delete flowcells')

        flowcell_list = Flowcell.objects.all()

        for flowcell in flowcell_list:

            flowcell_job_list = JobMaster.objects.filter(flowcell=flowcell)

            for flowcell_job in flowcell_job_list:

                if flowcell_job.job_type.name == "Delete_Flowcell":

                    print('Deleting flowcell {}'.format(flowcell.id))

                    flowcell.delete()
