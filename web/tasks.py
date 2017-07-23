from __future__ import absolute_import, unicode_literals
from celery import task

@task()
def task_number_one():
    # Do something...
    print ("Here I am running task_number_one")
    print ("get in", counter)