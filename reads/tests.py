import datetime

import pytz
from django.contrib.auth.models import User
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.test import APITestCase

from .models import Barcode, FastqRead, FastqReadType, Run, Flowcell


class MinionRunTestCase(APITestCase):

    def setUp(self):

        user_john = User.objects.create_user('john', 'john@test.com', 'passasdf')
        user_alex = User.objects.create_user('alex', 'alex@test.com', 'passasdf')

        flowcell1 = Flowcell.objects.create(
            name='FL0001',
            owner=user_john
        )

        flowcell2 = Flowcell.objects.create(
            name='FL0002',
            owner=user_alex
        )

        run1 = Run.objects.create(name='20170516_1630_john',
                                  runid='hj78yy9o-217e-4335-9451-66a7288a9dd5',
                                  is_barcoded=False,
                                  flowcell=flowcell1,
                                  owner=user_john)

        Run.objects.create(name='20170517_1630_john',
                           runid='hj78yy9o-217e-4335-9451-66a7288a9aa6',
                           is_barcoded=True,
                           flowcell=flowcell1,
                           owner=user_john)

        Run.objects.create(name='20170518_1630_alex',
                           runid='hj78yy9o-217e-4335-9451-66a7288a9aa6',
                           is_barcoded=True,
                           flowcell=flowcell2,
                           owner=user_alex)

        barcode1 = Barcode.objects.create(name='Barcode 1', run=run1)

        FastqRead.objects.create(run=Run.objects.get(name='20170516_1630_john'),
                                 read_id='d43c9b02-5d9f-4d37-9ed4-b718371d1f86',
                                 read=6632,
                                 channel=231,
                                 barcode=barcode1,
                                 sequence='ACATTTGCCG',
                                 quality='$#%*&)''%#',
                                 is_pass=True,
                                 type=FastqReadType.objects.get(name='Template'),
                                 start_time=datetime.datetime(2017, 5, 3, 16, 28, 5, 0, pytz.UTC))

        FastqRead.objects.create(run=Run.objects.get(name='20170516_1630_john'),
                                 read_id='2e48e4a3-aaf7-4311-91fe-4bb2054e3bba',
                                 read=6316,
                                 channel=172,
                                 barcode=barcode1,
                                 sequence='ACATTTGCCG',
                                 quality='$#%*&)''%#',
                                 is_pass=True,
                                 type=FastqReadType.objects.get(name='Template'),
                                 start_time=datetime.datetime(2017, 5, 4, 16, 28, 5, 0, pytz.UTC))

        FastqRead.objects.create(run=Run.objects.get(name='20170518_1630_alex'),
                                 read_id='2e48e4a3-aaf7-4311-91fe-4bb2054e3bba',
                                 read=6316,
                                 channel=172,
                                 sequence='ACATTTGCCG',
                                 quality='$#%*&)''%#',
                                 is_pass=True,
                                 type=FastqReadType.objects.get(name='Template'),
                                 start_time=datetime.datetime(2017, 5, 5, 16, 28, 5, 0, pytz.UTC))

        FastqRead.objects.create(run=Run.objects.get(name='20170518_1630_alex'),
                                 read_id='2e48e4a3-aaf7-4311-91fe-4bb2054e3bba',
                                 read=6316,
                                 channel=172,
                                 sequence='ACATTTGCCG',
                                 quality='$#%*&)''%#',
                                 is_pass=True,
                                 type=FastqReadType.objects.get(name='Template'),
                                 start_time=datetime.datetime(2017, 5, 6, 16, 28, 5, 0, pytz.UTC))

    def test_minionruns_have_barcode(self):
        minionrun_no_barcoded = Run.objects.get(run_name='20170516_1630_john')
        minionrun_barcoded = Run.objects.get(run_name='20170517_1630_john')

        self.assertEqual(minionrun_barcoded.is_barcoded, True)
        self.assertEqual(minionrun_no_barcoded.is_barcoded, False)

    def test_minionruns_get_method_with_wrong_auth(self):
        response = self.client.get('/api/v1/runs/')

        self.client.credentials(HTTP_AUTHORIZATION='Token a4ef48961c860e72a68534e29dc4ec8f2aea344f')

        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)
