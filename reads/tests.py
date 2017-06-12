import datetime
import pytz

from django.contrib.auth.models import User
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.test import APITestCase

from .models import MinIONRun, FastqReadType, FastqRead


class MinionRunTestCase(APITestCase):
    def setUp(self):
        user_john = User.objects.create_user('john', 'john@test.com', 'passasdf')
        user_alex = User.objects.create_user('alex', 'alex@test.com', 'passasdf')

        MinIONRun.objects.create(run_name='20170516_1630_john',
                                 run_id='hj78yy9o-217e-4335-9451-66a7288a9dd5',
                                 is_barcoded=False,
                                 owner=user_john)

        MinIONRun.objects.create(run_name='20170517_1630_john',
                                 run_id='hj78yy9o-217e-4335-9451-66a7288a9aa6',
                                 is_barcoded=True,
                                 owner=user_john)

        MinIONRun.objects.create(run_name='20170518_1630_alex',
                                 run_id='hj78yy9o-217e-4335-9451-66a7288a9aa6',
                                 is_barcoded=True,
                                 owner=user_alex)

        FastqReadType.objects.create(name='Template')
        FastqReadType.objects.create(name='Complement')
        FastqReadType.objects.create(name='2d')

        FastqRead.objects.create(run_id=MinIONRun.objects.get(run_name='20170516_1630_john'),
                                 read_id='d43c9b02-5d9f-4d37-9ed4-b718371d1f86',
                                 read=6632,
                                 channel=231,
                                 barcode='barcode1',
                                 sequence='ACATTTGCCG',
                                 quality='$#%*&)''%#',
                                 is_pass=True,
                                 type=FastqReadType.objects.get(name='Template'),
                                 start_time=datetime.datetime(2017, 5, 3, 16, 28, 5, 0, pytz.UTC))

        FastqRead.objects.create(run_id=MinIONRun.objects.get(run_name='20170516_1630_john'),
                                 read_id='2e48e4a3-aaf7-4311-91fe-4bb2054e3bba',
                                 read=6316,
                                 channel=172,
                                 barcode='barcode1',
                                 sequence='ACATTTGCCG',
                                 quality='$#%*&)''%#',
                                 is_pass=True,
                                 type=FastqReadType.objects.get(name='Template'),
                                 start_time=datetime.datetime(2017, 5, 4, 16, 28, 5, 0, pytz.UTC))

        FastqRead.objects.create(run_id=MinIONRun.objects.get(run_name='20170518_1630_alex'),
                                 read_id='2e48e4a3-aaf7-4311-91fe-4bb2054e3bba',
                                 read=6316,
                                 channel=172,
                                 barcode='barcode1',
                                 sequence='ACATTTGCCG',
                                 quality='$#%*&)''%#',
                                 is_pass=True,
                                 type=FastqReadType.objects.get(name='Template'),
                                 start_time=datetime.datetime(2017, 5, 5, 16, 28, 5, 0, pytz.UTC))

        FastqRead.objects.create(run_id=MinIONRun.objects.get(run_name='20170518_1630_alex'),
                                 read_id='2e48e4a3-aaf7-4311-91fe-4bb2054e3bba',
                                 read=6316,
                                 channel=172,
                                 barcode='barcode1',
                                 sequence='ACATTTGCCG',
                                 quality='$#%*&)''%#',
                                 is_pass=True,
                                 type=FastqReadType.objects.get(name='Template'),
                                 start_time=datetime.datetime(2017, 5, 6, 16, 28, 5, 0, pytz.UTC))

    def test_minionruns_have_barcode(self):
        minionrun_no_barcoded = MinIONRun.objects.get(run_name='20170516_1630_john')
        minionrun_barcoded = MinIONRun.objects.get(run_name='20170517_1630_john')

        self.assertEqual(minionrun_barcoded.is_barcoded, True)
        self.assertEqual(minionrun_no_barcoded.is_barcoded, False)

    def test_minionruns_get_method_without_auth(self):
        response = self.client.get('/api/v1/runs/')

        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_minionruns_get_method_with_wrong_auth(self):
        response = self.client.get('/api/v1/runs/')

        self.client.credentials(HTTP_AUTHORIZATION='Token a4ef48961c860e72a68534e29dc4ec8f2aea344f')

        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    """
    def test_minionruns_get_method_with_correct_auth(self):
        token = Token.objects.get(user__username='john')

        self.client.credentials(HTTP_AUTHORIZATION='Token ' + token.key)

        response = self.client.get('/api/v1/runs/')

        print(response.data)

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data,
                         [
                             {
                                 'url': 'http://192.168.100.10:8000/api/v1/runs/1/',
                                 'run_name': '20170516_1630_john',
                                 'run_id': 'hj78yy9o-217e-4335-9451-66a7288a9dd5',
                                 'is_barcoded': False
                             },
                             {
                                 'id': 'http://192.168.100.10:8000/api/v1/runs/2/',
                                 'run_name': '20170517_1630_john',
                                 'run_id': 'hj78yy9o-217e-4335-9451-66a7288a9aa6',
                                 'is_barcoded': True
                             }
                         ])

        token = Token.objects.get(user__username='alex')

        self.client.credentials(HTTP_AUTHORIZATION='Token ' + token.key)

        response = self.client.get('/api/v1/runs/')

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data,
                         [
                             {
                                 'id': 3,
                                 'run_name': '20170518_1630_alex',
                                 'run_id': 'hj78yy9o-217e-4335-9451-66a7288a9aa6',
                                 'is_barcoded': True
                             }
                         ])


    # this test doesn't work because the minionrun ids changed from 1 and 2 to 13 and 14.
    def test_runs_get_method_with_correct_owner(self):
        token = Token.objects.get(user__username='john')

        self.client.credentials(HTTP_AUTHORIZATION='Token ' + token.key)

        response = self.client.get('/api/v1/runs/')

        print(response.data)

        self.assertEqual(response.status_code, status.HTTP_200_OK)

        self.assertEqual(response.data,
                         [
                             {
                                 'id': 1,
                                 'run_name': '20170516_1630_john',
                                 'run_id': 'hj78yy9o-217e-4335-9451-66a7288a9dd5',
                                 'is_barcoded': False
                             }
                         ])
    """