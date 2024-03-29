from __future__ import absolute_import, unicode_literals

import datetime
import logging
import os
import pickle
import uuid
from pathlib import Path

from celery import Celery, Task
from celery.worker.request import Request

from minotourapp.settings import BASE_DIR

logger = logging.getLogger("my.package")


class MyRequest(Request):
    """A minimal custom request to log failures and hard time limits."""

    def on_timeout(self, soft, timeout):
        super(MyRequest, self).on_timeout(soft, timeout)
        error_log_file = f"{BASE_DIR}/celery_error.log"
        with open(error_log_file, "a") as fh:
            fh.write(
                f"{self.task.name} failed at {datetime.datetime.now()}"
                f" processing:\n {str(self.info())}"
            )

    def on_failure(self, exc_info, send_failed_event=True, return_ok=False):
        super(MyRequest, self).on_failure(
            exc_info, send_failed_event=send_failed_event, return_ok=return_ok
        )
        error_log_file = f"{BASE_DIR}/celery_error.log"
        uuid_str = uuid.uuid4()
        Path(f"{BASE_DIR}/pickled_args/").mkdir(exist_ok=True, parents=True)
        pickle_file = f"{BASE_DIR}/pickled_args/{uuid_str}.pickle"
        with open(error_log_file, "a") as fh:
            fh.write(
                f"{self.task.name} failed at {datetime.datetime.now()} with {str(exc_info.traceback)}"
                f" processing:\n args stored in pickle at {uuid_str}.pickle\n"
            )
        with open(pickle_file, "wb") as fh:
            pickle.dump(self._args, fh, protocol=pickle.HIGHEST_PROTOCOL)


class MyTask(Task):
    Request = MyRequest


# set the default Django settings module for the 'celery' program.
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "minotourapp.settings")

app = Celery("minotourapp")


# Using a string here means the worker don't have to serialize
# the configuration object to child processes.
# - namespace='CELERY' means all celery-related configuration keys
#   should have a `CELERY_` prefix.
app.config_from_object("django.conf:settings", namespace="CELERY")

app.autodiscover_tasks()
