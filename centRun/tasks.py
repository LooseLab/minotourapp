# create your tasks here, should return asynchronously

from __future__ import absolute_import, unicode_literals
from centRun.RunCentrifuge import Centrifuger
import sys
from django.http import HttpResponse
# from config.celery import app
from celery import shared_task


@shared_task
def start_centrifuge(foreign_key):
    print("starting centrifuge task")
    print(f"identity is {foreign_key}")
    """
    starts the centrifuge instance as a celery task which should return asynchronously
    :param request:
    :return:
    """
    try:
        print("trying centrifuge task")
        c = Centrifuger(foreign_key)
        c.run_centrifuge()
        print(sys.stdout)

    except:
        e = sys.exc_info()
        print(e)



