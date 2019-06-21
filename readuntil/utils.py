"""
Useful functions for this app
"""
from pathlib import Path


def make_results_directory(flowcell_id, task_id):
    """
    Make directory to store results of this task
    :param flowcell_id: The ID of the flowcell that produced these fastqs
    :param task_id: The ID of the Jobmaster that corresponds to this EB task
    :return: The path into this directory
    """
    results_dir = Path(f"readuntil/Temp_results/EB{flowcell_id}_{task_id}")
    if not results_dir.exists():
        Path.mkdir(results_dir)
    return results_dir