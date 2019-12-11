import logging

from django.test import TestCase

from reads.models import Flowcell, Run, FastqRead, JobType, JobMaster
from web.tasks import update_flowcell_details

logger = logging.getLogger(__name__)
logging.disable(logging.NOTSET)
logger.setLevel(logging.DEBUG)


class UpdateFlowcellDetails(TestCase):

    fixtures = ['fixtures/auxiliary_data.json', 'fixtures/reference.json', 'fixtures/yeast.json']

    def test_update_flowcell_details(self):

        print('Test 122')

        flowcell_list = Flowcell.objects.all()
        flowcell = flowcell_list[0]

        run_list = Run.objects.filter(flowcell=flowcell)

        for run in run_list:
            reads_list = FastqRead.objects.filter(run=run)

            print('Run {} - {} - reads {}'.format(run.id, run, len(reads_list)))

        job_type_readuntil = JobType.objects.get(pk=15)
        job_master = JobMaster.objects.create(
            flowcell=flowcell,
            job_type=job_type_readuntil
        )

        update_flowcell_details(job_master.id)

        print('End')

        # self.assertEqual(line_tuple.read_id, '3d8564a8-653e-4dbe-b301-3a2cea209bf8')
        # self.assertEqual(line_tuple.chromosome, 'NC_001140')

