from django.core.management import BaseCommand, CommandError

from reads.models import Flowcell
from jobs.models import JobMaster
from web.tasks import update_flowcell_details, processreads, run_centrifuge


class Command(BaseCommand):

    help = 'Run task centrifuge'

    def handle(self, *args, **options):

        try:

            print('Running centrifuge task')

            flowcell_list = Flowcell.objects.all()

            for flowcell in flowcell_list:

                flowcell_job_list = JobMaster.objects.filter(flowcell=flowcell)

                for flowcell_job in flowcell_job_list:

                    if flowcell_job.job_type.name == "Metagenomics":

                        print("trying to run classification for flowcell {} {} {} ".format(
                            flowcell.id,
                            flowcell_job.id,
                            flowcell_job.last_read
                        ))

                        print("starting centrifuge task")

                        run_centrifuge(flowcell.id, flowcell_job.id)


        except Exception as e:

            raise CommandError(repr(e))