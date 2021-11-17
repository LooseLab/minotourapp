import re
from pathlib import Path
from typing import Union

import numpy as np
from pyfastx import Fastx

from alignment.models import PafSummaryCov
from minotourapp.utils import get_env_variable
from reads.models import Barcode


def _human_key(key):
    parts = re.split("(\d*\.\d+|\d+)", key)
    return tuple(
        (e.swapcase() if i % 2 == 0 else float(e)) for i, e in enumerate(parts)
    )


def read_fastq_to_dict(file_path):
    """
    File path
    Parameters
    ----------
    file_path: str
        The file path to the fastq file
    Returns
    -------
    dict
        Fastq dictionary
    """
    fx = Fastx(file_path)
    result_list = []
    for name, seq, qual, comment in fx:
        comment_dict = dict(map(lambda x: x.split("="), comment.split()))
        result_list.append(
            {
                "read_id": name,
                "barcode_name": comment_dict.get("barcode", ""),
                "sequence": seq,
                "run_id": comment_dict.get("runid"),
                "type__name": "T",
            }
        )
    return result_list


def get_alignment_result_dir(run_id: str, create: bool) -> Path:
    """
    Get the alignment folder for this run
    Parameters
    ----------
    run_id: str
        The run id of the data
    create: bool
        Create the directory if it doesn't exist
    Returns
    -------
    Path
        Path to the Artic results dir
    """
    alignment_results_path = (
        Path(get_env_variable("MT_ALIGNMENT_DATA_DIR")) / "alignment" / run_id
    )
    if create:
        alignment_results_path.mkdir(exist_ok=True, parents=True)
    if not alignment_results_path.exists():
        raise FileNotFoundError(f"{alignment_results_path} does not exist.")
    return alignment_results_path


def get_or_create_array(
    folder_dir: Path,
    contig_name: str,
    contig_length: int = None,
    barcode_name: Union[bool, str] = False,
    bin_width: int = 10,
    create: bool = False,
    is_cnv: bool = False
) -> Path:
    """
    Get the path to the numpy array containing the binned alignment
     counts, create it if it doesn't exist. Created as an array of numpy uint16,
     with length of contig // 10
    Parameters
    ----------
    folder_dir: Path
        The path to the folder containing these alignments
    contig_name: str
        The name of the contig for the alignment
    contig_length: int, default None
        The length of the contig, used if creating the array
    barcode_name: bool or str, default False
        The name of the barcode if barcoded run, else False
    bin_width: int
        The width of the bins to use
    create: bool
        Create the array if it doens't exists
    is_cnv: bool
        Return the coverage pile up bins
    Returns
    -------
    Path
    """
    # only need 1 dimension if is
    shape_me = 1 if is_cnv else 2
    barcode_dir = "no_barcode" if not barcode_name else barcode_name
    suffice = "cnv_bins" if is_cnv else "bins"
    array_path = Path(f"{folder_dir/barcode_dir/contig_name}/{contig_name}_{suffice}.npy")
    if not array_path.exists() and create:
        dtype_me_baby = np.uint16 if contig_length < 100_000_000 else np.uint8
        array_path.parent.mkdir(parents=True, exist_ok=True)
        array = np.zeros((shape_me, np.ceil(contig_length / bin_width).astype(int)), dtype=np.uint16)
        np.save(array_path, arr=array)
    return array_path


def create_paf_summary_cov(job_master, read_dict: dict, reference, contig):
    """
    Create a paf summary cov
    Parameters
    ----------
    job_master: read.models.JobMaster
        The job master for this alignment task
    read_dict: dict
        Dictionary containing all the read_info we need
    reference: reference.models.ReferenceInfo
        The reference being mapped against
    contig: reference.models.ReferenceLine
        The contig that we are mapping against

    Returns
    -------
    alignment.models.PafSummaryCov
    """
    barcode = Barcode.objects.get(pk=read_dict["barcode"])
    return PafSummaryCov.objects.get_or_create(
        job_master=job_master,
        barcode_id=read_dict["barcode"],
        barcode_name=barcode.name,
        chromosome=contig,
        reference_line_length=contig.chromosome_length,
        read_type_id=read_dict["type"],
        chromosome_pk=contig.id,
        chromosome_name=contig.line_name,
        reference_pk=reference.id,
        reference_name=reference.name,
    )[0]
