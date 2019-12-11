import logging
import os

from django.core.management import BaseCommand, CommandError

from reads.models import Flowcell, FastqRead, JobMaster

log_folder = os.environ.get('MT_LOG_FOLDER')

logging.basicConfig(
    format='%(asctime)s %(module)s:%(levelname)s:%(thread)d:%(message)s',
    filename=os.path.join(log_folder, 'run_task_delete_flowcells.log'),
    level=os.environ.get('LOGLEVEL', 'INFO')
)

logger = logging.getLogger(__name__)


class Command(BaseCommand):

    help = 'Run task delete flowcells'

    def handle(self, *args, **options):

        logger.info('Running delete flowcells task')

        try:

            flowcell_list = Flowcell.objects.order_by('id').values('id', 'name')

            for flowcell in flowcell_list:

                flowcell_job_list = JobMaster.objects.filter(flowcell__id=flowcell['id']).values('id', 'job_type__name')

                for flowcell_job in flowcell_job_list:

                    if flowcell_job['job_type__name'] == "Delete_Flowcell":

                        logger.info('Flowcell id: {} - Deleting flowcell {}'.format(flowcell['id'], flowcell['name']))

                        first_fastqread = FastqRead.objects.filter(run__flowcell__id=flowcell['id']).order_by('id').first()
                        last_fastqread = FastqRead.objects.filter(run__flowcell__id=flowcell['id']).order_by('id').last()

                        if first_fastqread and last_fastqread:

                            first_read = first_fastqread.id
                            last_read = last_fastqread.id

                            logger.info('Flowcell id: {} - First and last fastqread id are {} {}'.format(flowcell['id'], first_read, last_read))

                            counter = first_read

                            while counter <= last_read:

                                counter = counter + 50

                                logger.info('Flowcell id: {} - Deleting records with id < {}'.format(flowcell['id'], counter))

                                affected = FastqRead.objects.filter(run__flowcell__id=flowcell['id']).filter(id__lt=counter).delete()

                                logger.info('Flowcell id: {} - Deleted {} fastqread records'.format(flowcell['id'], affected))

                        affected = Flowcell.objects.get(pk=flowcell['id']).delete()

                        logger.info('Flowcell id: {} - Deleted {} records'.format(flowcell['id'], affected))

        except Exception as e:

            logger.exception('Failed running task.')
            raise CommandError(repr(e))

        finally:

            logger.info('Finishing delete flowcells task')
