from typing import List

from minotourapp.celery import app


@app.task
def task_cnv(mappings: List[List], bin_width: int = 10) -> None:
    """
    Bin the mapping starts into bins of `bin_width`
    Parameters
    ----------
    mappings: list
        The mappings, as produced from the Mapper function
    bin_width: int default 10
        The width of the bin to bin into
    Returns
    -------
    None
    """
    read_ids = {""}
    (mapping for mapping in mappings)
