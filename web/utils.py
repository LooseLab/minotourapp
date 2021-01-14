import datetime
import os

from django.contrib.auth.models import User

from communication.models import Message
from minknow_data.models import MinionRunInfo, Flowcell
from minotourapp.utils import get_env_variable


def file_folder_exists(filepath):
    """
    Check if a provided file or folder exists.

    Parameters
    ----------
    filepath: str
        Absolute path to file or folder to be tested.

    Returns
    -------
    bool

    """
    return os.path.exists(filepath)


def is_executable(filepath):
    """
    Check whether the provided file path is executable.
    Parameters
    ----------
    filepath: str

    Returns
    -------
    bool

    """
    return os.access(filepath, os.X_OK)


def validate_envs():
    """
    Validate environment variables for obvious errors.
    Returns
    -------
    None
    TODO finish this
    """
    # Artic checks
    # if " " in get_env_variable("MT_ARTIC_RESULTS_DIR"):
    #     raise ValueError("Please remove whitespace from the Artic results directory variable name.")
    # if not file_folder_exists(get_env_variable("MT_ARTIC_RESULTS_DIR")):
    #     raise ValueError("The Artic results directory specified in the environment variables file does not exist.")
    #
    # # Minimap2 checks
    # if not is_executable(get_env_variable("MT_MINIMAP2")):
    #     raise ValueError("The Minimap2 executable is either not there, or not executable. Please double check.")
    pass


def update_last_activity_time(flowcell):
    """
    Update the last activity date on the flowcell
    :param flowcell: The flowcell to be updated
    :return:
    """
    flowcell.last_activity_date = datetime.datetime.now(datetime.timezone.utc)
    flowcell.save()


def get_run_details(run_id):
    # TODO i hate this code, let's rewrite this code and comment it - SERIOUSLY HORRENDOUS
    flowcell = Flowcell.objects.get(pk=run_id)

    result = {}

    for run in flowcell.runs.all():

        element = {
            "id": None,
            "runid": None,
            "run_start_time": None,
            "first_read": None,
            "last_read": None,
            "minknow_computer_name": None,
            "minion_id": None,
            "asic_id": None,
            "sequencing_kit": None,
            "purpose": None,
            "minknow_version": None,
            "flowcell_type": None,
            "flowcell_id": None,
            "sample_name": None,
            "experiment_name": None,
            "read_count": None,
            "total_read_length": None,
            "max_read_length": None,
            "min_read_length": None,
            "avg_read_length": None,
            "first_read_start_time": None,
            "last_read_start_time": None,
        }

        minion_run_status_list = MinionRunInfo.objects.filter(run_id=run).order_by(
            "minKNOW_start_time"
        )

        if len(minion_run_status_list) > 0:

            minion_run_status = minion_run_status_list[0]

            element["id"] = minion_run_status.run.id
            element["runid"] = minion_run_status.run.runid
            element["run_start_time"] = minion_run_status.minKNOW_start_time
            element["minknow_computer_name"] = minion_run_status.minKNOW_computer
            element["minion_id"] = minion_run_status.minion.minION_name
            element["asic_id"] = minion_run_status.minKNOW_asic_id
            element["sequencing_kit"] = minion_run_status.sequencing_kit
            element["purpose"] = minion_run_status.minKNOW_exp_script_purpose
            element["minknow_version"] = minion_run_status.minKNOW_version
            element["flowcell_type"] = minion_run_status.flowcell_type
            element["flowcell_id"] = minion_run_status.minKNOW_flow_cell_id
            element["sample_name"] = minion_run_status.minKNOW_sample_name
            element["experiment_name"] = minion_run_status.experiment_id

            result[element["runid"]] = element

    result_basecalled_summary = []

    for run in flowcell.runs.all():

        if hasattr(run, "summary"):

            runid = run.summary.run.runid

            if runid in result.keys():

                result[runid]["read_count"] = run.summary.read_count
                result[runid]["total_read_length"] = run.summary.total_read_length
                result[runid]["max_read_length"] = run.summary.max_read_length
                result[runid]["min_read_length"] = run.summary.min_read_length
                result[runid]["avg_read_length"] = run.summary.avg_read_length
                result[runid][
                    "first_read_start_time"
                ] = run.summary.first_read_start_time
                result[runid]["last_read_start_time"] = run.summary.last_read_start_time

            else:

                element = {
                    "runid": run.summary.run.runid,
                    "run_start_time": None,
                    "first_read": None,
                    "last_read": None,
                    "minknow_computer_name": None,
                    "minion_id": None,
                    "asic_id": None,
                    "sequencing_kit": None,
                    "purpose": None,
                    "minknow_version": None,
                    "flowcell_type": None,
                    "flowcellId": None,
                    "sample_name": None,
                    "experiment_name": None,
                    "read_count": run.summary.read_count,
                    "total_read_length": run.summary.total_read_length,
                    "max_read_length": run.summary.max_read_length,
                    "min_read_length": run.summary.min_read_length,
                    "avg_read_length": run.summary.avg_read_length,
                    "first_read_start_time": run.summary.first_read_start_time,
                    "last_read_start_time": run.summary.last_read_start_time,
                    "id": run.id,
                }

                result[element["runid"]] = element

    return result.values()


def fun(self, exc, task_id, args, kwargs, einfo):
    """
    On error function for tasks - gonna set up so admin account can
    Parameters
    ----------
    self
    exc: Exception raised by the task
    task_id: Unique Id of the failed task
    args: original arguments of the faield task
    kwargs: keyword arguments of the failed task
    einfo: ExceptionInfo
        Info about the exception raised

    Returns
    -------

    """
    u = User.objects.filter(is_superuser=True)[0]
    if u and int(get_env_variable("MT_TWEET_CELERY_FAILURES")):
        error_string = einfo.__str__().split("\n")[-2]
        Message.objects.create(
            recipient=u,
            title=f"Celery task {self.name} failed at"
            f" {str(datetime.datetime.now())}. Final stack trace line: "
            f"{error_string} ",
        )
