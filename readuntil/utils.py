"""
Useful functions for this app
"""
from pathlib import Path
from django.conf import settings

def get_or_create_results_directory(flowcell_id, task_id):
    """
    Make directory to store results of this task
    :param flowcell_id: The ID of the flowcell that produced these fastqs
    :param task_id: The ID of the Jobmaster that corresponds to this EB task
    :return: The path into this directory
    """
    current_working_directory = getattr(settings, "REFERENCE_LOCATION", None)
    results_dir = Path(f"{current_working_directory}/readuntil/Temp_results/EB{flowcell_id}_{task_id}")
    if not results_dir.exists():
        Path.mkdir(results_dir)
    return results_dir


def check_results_exist(flowcell_id, task_id, chromosome):
    """
    Check if we have results files present.
    :param flowcell_id: The primary of key of the flowcell database entry
    :type flowcell_id: int
    :param task_id: The primary key of the task database entry
    :type task_id: int
    :param chromosome: The name of a chromosome we have performed the readUntil task on.
    :type chromosome: str
    :return: boolean confirming if we have a results file or not
    """

    results_dir = get_or_create_results_directory(flowcell_id, task_id)

    coverage_path = Path(
        f"{results_dir}/coverage_{chromosome}_{flowcell_id}_{task_id}.dat"
    )

    exists = coverage_path.exists()

    return exists
