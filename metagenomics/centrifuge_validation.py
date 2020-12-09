import gzip
import subprocess
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import pyfastx
from celery.utils.log import get_task_logger
from django.conf import settings
from ete3 import NCBITaxa

from alignment.tasks_alignment import run_minimap2_alignment
from metagenomics.models import (
    MappingTarget,
    MappingResult,
    CentrifugeOutput,
)
from metagenomics.utils import falls_in_region
from minotourapp.utils import get_env_variable
from reads.models import JobMaster
from reference.models import ReferenceInfo, ReferenceLine
from reference.utils import create_minimap2_index

logger = get_task_logger(__name__)


def separate_target_cent_output(df, task, fasta_df_barcode):
    """
    Separate out any reads from the dataframe that classified as targets
    Parameters
    ----------
    df: pandas.core.frame.DataFrame
         The metagenomics output dataframe
    task: reads.models.JobMaster
        The task ORM that is representing this metagenomics task
    Returns
    -------

    """

    targets = (
        MappingTarget.objects.filter(target_set=task.target_set)
        .values_list("tax_id", flat=True)
        .distinct()
    )
    logger.info("Flowcell id: {} - Targets are {}".format(task.flowcell.id, targets))
    # Get the target reads rows in a separate dataframe
    targets_cent_out_df = df[df["tax_id"].isin(targets)]
    target_species_df = pd.merge(
        fasta_df_barcode, targets_cent_out_df, on="read_id", how="inner"
    )
    target_species_df.drop(columns=["barcode_name_y"], inplace=True)
    target_species_df.rename(columns={"barcode_name_x": "barcode_name"}, inplace=True)
    # Get the results into  a dataframe
    return target_species_df


def create_mapping_result_objects(barcodes, task):
    """
    Create the mapping results objects for the all reads barcode and any barcode that is in this iteration
    Parameters
    ----------
    barcodes: list of str
        The barcodes that we have identified reads for, plus ALl reads
    task: reads.models.JobMaster
        The Django ORM object for this task

    Returns
    -------
    pandas.core.frame.DataFrame
        A dataframe of mapping target information

    """
    target_regions_df = pd.DataFrame.from_records(
        MappingTarget.objects.filter(target_set=task.target_set).values()
    )
    # If there are no target regions found for the set name
    if target_regions_df.empty:
        logger.info(f"No target set pk {task.target_set}")
        return
    # Get np array of the target species
    target_species = target_regions_df["species"].unique()
    # Get np array of their tax_ids
    target_tax_id = target_regions_df["tax_id"].unique()
    # Combine into tuple, one for each
    target_tuple = list(zip(target_species, target_tax_id))
    # create one MappingResults entry for each target species
    for target in target_tuple:
        for barcode in barcodes:
            obj, created = MappingResult.objects.get_or_create(
                task=task,
                species=target[0],
                tax_id=target[1],
                barcode_name=barcode,
                defaults={"red_reads": 0, "num_mapped": 0},
            )
            if created:
                logger.info(
                    f"Flowcell id: {task.flowcell.id} - Result Mapping object created for {obj.species}, barcode {barcode}"
                )
    return target_regions_df


def create_concat_reference_info(task, concat_file_path):
    """
    Create a reference info database entry for the concatenated reference
    Parameters
    ----------
    task: reads.models.JobMaster
        Django ORM object for the JobMaster
    concat_file_path: pathlib.PosixPath
        The file path to the concatenated reference file.
    Returns
    -------
    None

    """
    ref_info, created = ReferenceInfo.objects.get_or_create(
        name=task.target_set,
        file_location=str(concat_file_path.resolve()),
        file_name=concat_file_path.stem,
    )
    handle = (
        pyfastx.Fasta
        if set(concat_file_path.suffixes).intersection(
            {".fna", ".fa", ".fsa", ".fasta"}
        )
        else pyfastx.Fastq
    )
    fa = handle(concat_file_path.as_posix())
    for contig in fa:
        ReferenceLine.objects.create(
            reference=ref_info, line_name=contig.name, chromosome_length=len(contig)
        )
    minimap2_index_path = create_minimap2_index(ref_info, concat_file_path)
    ref_info.minimap2_index_file_location = minimap2_index_path
    ref_info.save()


def fetch_concat_reference_path(task):
    """
    Fetch the path to the concat references of all members of a target set.
     Create it if doesn't exist.
    Parameters
    ----------
    task: reads.models.JobMaster
        The Django ORM object for this task

    Returns
    -------
    pathlib.PurePath
        Path to the concatenated references of a target set

    """
    concat_reference_path = Path(
        f"{settings.MEDIA_ROOT}/reference_files/{task.target_set}.fna.gz"
    ).resolve()
    if not concat_reference_path.exists():
        logger.info(
            f"Concatenated reference for target set {task.target_set} does not exist."
            f" Creating at location {concat_reference_path.as_posix()}"
        )
        reference_pks = set(
            MappingTarget.objects.filter(target_set=task.target_set).values_list(
                "reference_id", flat=True
            )
        )
        with gzip.open(concat_reference_path, "wb") as fh:
            for pk in reference_pks:
                ref_file = ReferenceInfo.objects.get(pk=pk)
                handle = (
                    gzip.open
                    if Path(ref_file.file_location.path).suffix == ".gz"
                    else open
                )
                fh.write(
                    handle(ref_file.file_location.path, "rt").read().strip().encode()
                )
                fh.write(b"\n")
        create_concat_reference_info(task, concat_reference_path)
    return concat_reference_path


def start_metagenomics_mapping_task(flowcell_id, target_set, target_reads_df):
    """
    Start a job type of metagenomics mapping that calls the alignment pipeline on this data
    Parameters
    ----------
    flowcell_id: int
        The primary key of the flowcell record in the database
    target_set: str
        The target set that we are working on
    target_reads_df: pd.core.frame.DataFrame
        Target reads dataframe of sequence read_id etc.
    Returns
    -------
    None
    """
    # TODO need a way of rerunning metagenomics analyses
    reference_info = ReferenceInfo.objects.get(name=target_set)
    jm, created = JobMaster.objects.get_or_create(
        job_type_id=19,
        flowcell_id=flowcell_id,
        target_set=target_set,
        reference_id=reference_info.id,
    )
    run_minimap2_alignment.apply_async(
        args=(jm.id, target_reads_df.to_dict(orient="records")), queue="minimap"
    )


def map_target_reads(task, path_to_reference, target_df, to_save_df, target_region_df):
    """
    Map the reads
    Parameters
    ----------
    task: reads.models.JobMaster
        The django ORM object of this task
    path_to_reference: pathlib.PosixPath
        The path to the concatenated, Gzipped reference file for all references in target set
    target_df: pd.core.frame.DataFrame
        Target reads dataframe, containing read sequence
    to_save_df: pd.core.frame.DataFrame
        The finalised metagenomics output data, with num_matches, lineages etc.
    target_region_df: pd.core.frame.DataFrame
        The regions defined in the GFF file that contain virulence areas

    Returns
    -------
    pd.core.frame.DataFrame
        Dataframe of minimap2 output from mapped target reads and num_mapped
    """
    minimap2_executable_path = get_env_variable("MT_MINIMAP2")
    cmd = f"{minimap2_executable_path} -x map-ont {path_to_reference} -"
    # target_df = pd.merge(target_df, to_save_df, on="tax_id")
    target_df["unique"] = np.where(target_df["numMatches"] == 1, 1, 0)
    gb = target_df.groupby(["tax_id", "barcode_name"])
    target_df.set_index(["tax_id", "barcode_name"], inplace=True)
    target_df["num_matches"] = gb.size()
    target_df["sum_unique"] = gb["unique"].sum()
    target_df.reset_index(inplace=True)
    taxid_list = np.unique(target_df["tax_id"].values)
    # use the ncbi thing to get the species names for the tax ids in the tax ids list
    ncbi = NCBITaxa()
    taxid_2_name = ncbi.get_taxid_translator(taxid_list)
    target_df["name"] = target_df["tax_id"].map(taxid_2_name)
    fasta_sequence_to_map = "\n".join(
        (">" + target_df["read_id"] + "\n" + target_df["sequence"]).values.tolist()
    )
    start_metagenomics_mapping_task(task.flowcell.id, task.target_set, target_df)
    # TODO merge mapping task into this function rather than above seperate mapping task creation
    process = subprocess.Popen(
        cmd.split(),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    out, err = process.communicate(input=fasta_sequence_to_map.encode())
    if not out:
        logger.info(f"No reads mapped.")
        return pd.DataFrame()
    target_set_plasmids_ref_pks = MappingTarget.objects.filter(
        target_set=task.target_set
    ).values_list("reference_id", flat=True)
    reference_contig_names = (
        ReferenceInfo.objects.filter(pk__in=target_set_plasmids_ref_pks)
        .values_list("reference_lines__line_name", "name")
        .distinct()
    )
    contig_to_reference_dict = {
        contig_name: ref_species for contig_name, ref_species in reference_contig_names
    }
    map_out_df = pd.read_csv(StringIO(out.decode()), sep="\t", header=None)
    map_out_df.rename(
        columns={
            0: "read_id",
            1: "query_seq_len",
            2: "query_start",
            3: "query_end",
            4: "rel_strand",
            5: "target_seq_name",
            6: "target_seq_length",
            7: "target_start",
            8: "target_end",
            9: "num_residue_matches",
            10: "alignment_block_length",
            11: "mapping_qual",
        },
        inplace=True,
    )
    map_out_df["name"] = map_out_df["target_seq_name"].map(contig_to_reference_dict)
    map_out_df["num_residue_matches"] = map_out_df["num_residue_matches"].astype(
        np.int64
    )
    # Filter out low quality matches or the
    map_out_df = map_out_df.query("mapping_qual >= 40 & num_residue_matches >= 200")
    if map_out_df.empty:
        logger.info("Insufficient quality mappings.")
        return pd.DataFrame()
    # See whether the start or end of a mapping falls into the region
    map_out_df["read_is_red"] = 0
    map_out_df["read_is_red"] += target_region_df.apply(
        falls_in_region, args=(map_out_df,), axis=1
    )
    map_out_df["read_is_red"] = np.where(map_out_df["read_is_red"], 1, 0)
    map_out_df["name"] = map_out_df["name"].str.replace("_", " ")
    map_out_df = pd.merge(map_out_df, target_df, how="left", on=["read_id", "name"])
    map_out_df["barcode_name"] = map_out_df["read_id"].map(
        target_df.set_index("read_id")["barcode_name"].loc[
            ~target_df.set_index("read_id")["barcode_name"].index.duplicated()
        ]
    )
    map_out_df = map_out_df.fillna(0)
    map_out_df["barcode_name"] = np.where(
        map_out_df["barcode_name"] == "No_barcode",
        "All reads",
        map_out_df["barcode_name"],
    )
    gb = map_out_df.groupby(["barcode_name", "name"])
    map_out_df.set_index(["barcode_name", "name"], inplace=True)
    map_out_df["num_mapped"] = gb.size()
    map_out_df["num_red_reads"] = gb["read_is_red"].sum()
    map_out_df = map_out_df.loc[~map_out_df.index.duplicated()]
    map_out_df.reset_index(inplace=True)
    return map_out_df

    # target_df[]


def save_mapping_results(row, task):
    """
    Save the mapping results
    Returns
    -------

    """
    species_info = CentrifugeOutput.objects.filter(
        task=task, species=row["name"], barcode_name=row["barcode_name"]
    ).last()
    num_matches, sum_unique = (
        (species_info.num_matches, species_info.sum_unique)
        if species_info
        else (row["num_matches"], row["sum_unique"])
    )
    species = " ".join(row["name"].split(" ")[:2])
    map_result = MappingResult.objects.get(
        task=task, species=species, barcode_name=row["barcode_name"]
    )
    map_result.num_mapped += row["num_mapped"]
    map_result.num_matches += num_matches
    map_result.sum_unique += sum_unique
    map_result.red_reads += row["num_red_reads"]
    map_result.save()
