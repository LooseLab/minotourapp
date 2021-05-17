from __future__ import absolute_import, unicode_literals

import gzip
import os
import pickle
import subprocess
from collections import defaultdict
from copy import deepcopy
from datetime import datetime, timezone, timedelta
from io import StringIO
from shutil import copy, rmtree

from celery.schedules import crontab
from celery.utils.log import get_task_logger
from django.conf import settings
from git import Repo

from alignment.models import PafSummaryCov
from artic.models import ArticBarcodeMetadata, ArticFireConditions
from artic.task_write_out_artic_metrics import (
    write_out_artic_metrics,
    get_or_create_metrics_df,
)
from artic.utils import (
    get_amplicon_band_data,
    make_results_directory_artic,
    get_amplicon_stats,
    predict_barcode_will_finish,
)
from communication.models import Message
from minknow_data.models import Flowcell
from minotourapp.celery import MyTask
from minotourapp.celery import app
from minotourapp.settings import BASE_DIR, STATIC_ROOT
from minotourapp.utils import get_env_variable
from reads.models import (
    JobMaster,
    JobType,
    Barcode,
    FlowcellSummaryBarcode,
    FastqRead,
)
from readuntil.functions_EB import *

logger = get_task_logger(__name__)



@app.on_after_finalize.connect
def setup_periodic_tasks(sender, **kwargs):
    #check for update of VoCs every 24 hours.
    sender.add_periodic_task(
        crontab(hour=0, minute=0),
        #test.s('Happy Mondays!'),
        Update_VoCs.s(),
    )

@app.task
def Update_VoCs():
    MT_VoC_PATH = get_env_variable("MT_VoC_PATH")
    if Path(f"{MT_VoC_PATH}").exists():
        print("VoC Path Found")
        ##Check if
        # cloned_repo = Repo.clone(os.path.join("https://github.com/phe-genomics/variant_definitions", Path(f"{MT_VoC_PATH}")))
        if Path(f"{MT_VoC_PATH}/variant_definitions/").exists():
            # already cloned so....
            print("Updating path")
            repo = Repo(Path(f"{MT_VoC_PATH}/variant_definitions/"))
            print(repo.remotes.origin.pull())
            pass
        else:
            cloned_repo = Repo.clone_from("https://github.com/phe-genomics/variant_definitions",
                                          f"{MT_VoC_PATH}/variant_definitions/")


def check_afc_values_met(afc, coverage_array, num_amplicons, run_id, barcode_name, task):
    """
    Check if this artic fire condition has been met on this barcode
    Parameters
    ----------
    afc: artic.models.ArticFireConditions
        The model that we are checking
    coverage_array: np.ndarray
        Array with the mean coverage of each amplicon
    num_amplicons: int
        The number of amplicons
    run_id: int
        The primary key of this run
    barcode_name: str
        The name of this barcode as a string
    task: reads.models.JobMaster
        The job master for this artic task
    Returns
    -------
    bool
        Whether the Artic Firing Condition has been met

    """
    try:
        barcode = Barcode.objects.get(run_id=run_id, name=barcode_name)
        abm, created = ArticBarcodeMetadata.objects.get_or_create(flowcell=task.flowcell, job_master=task, barcode=barcode)
    except Barcode.DoesNotExist as e:
        logger.warning(f"Failed to update the Barcode metadata error: {e}")
    # check if the percent of amplicons median coverage is more than our preset and we haven't already fired
    # on this condition
    met_for_this_afc = (
        coverage_array[coverage_array > afc.x_coverage].size / num_amplicons * 100
    ) >= afc.percent_of_amplicons and afc not in abm.artic_fire_conditions.all()
    if met_for_this_afc:
        # we are firing on this condition so add it to the ABM so we don't fire on it again
        abm.artic_fire_conditions.add(afc)
    return met_for_this_afc


def get_amplicon_infos():
    """
    get amplicon band coordinates
    Returns
    -------
    (num_amplicons, amplicon_band_coords): (int, list)
        Tuple containing number of amplicons and the amplicon band coordinates

    """
    # todo hardcoded scheme atm, need to attach to JobMaster somehow
    scheme = "nCoV-2019"
    scheme_version = "V3"
    amplicon_band_coords, colours = get_amplicon_band_data(scheme, scheme_version)
    num_amplicons = len(amplicon_band_coords)
    return (num_amplicons, amplicon_band_coords)


# TODO parsing cigar lends itself to numba nicely
# TODO split artic task into a sexy chain, rather than blocking a worker for so long, start here
def clear_unused_artic_files(artic_results_path, sample_name, flowcell_id):
    """
    Clear the leftover files from the artic command from the directory, gzip remaining files
    Parameters
    ----------
    artic_results_path: str
        Path to directory that artic results were produced in.
    sample_name: str
        The sample name provided to the artic command. Assume the barcode name.
    flowcell_id: int
        The flowcell primary key for this artic task, so we can check that is has reads uploaded recently
    Returns
    -------
    None
    """
    flowcell = Flowcell.objects.get(pk=flowcell_id)
    logger.debug(
        f"Clearing non necessary results for {artic_results_path} for sample name {sample_name}"
    )
    files_to_keep_extra = [".fastq", ".sorted.bam", ".sorted.bam.bai"]
    files_to_keep = [
        ".pass.vcf",
        ".fail.vcf",
        "-barplot.png",
        "-boxplot.png",
        ".consensus.fasta",
        ".muscle.out.fasta",
        ".counts.dat",
        ".coverage.dat",
        "_ARTIC_medaka.csv"
    ]
    if int(
        get_env_variable("MT_DESTROY_ARTIC_EVIDENCE")
    ) and flowcell.last_activity_date < datetime.now(timezone.utc) - timedelta(
        hours=12
    ):
        logger.debug(
            f"Clearing sensitive files from {artic_results_path} for sample name {sample_name}"
        )
    else:
        logger.debug("Extending files to keep")
        files_to_keep.extend(files_to_keep_extra)
    files_to_keep_full = []
    files_to_keep_full.extend([f"{sample_name}{filey}" for filey in files_to_keep])
    files_to_keep_full.extend([f"{el}.gz" for el in files_to_keep_full])
    files_to_keep_full.append("lineage_report.csv")
    files_to_keep_full.append("lineage_report.csv.gz")
    files_to_keep_full.append("json_files/{sample_name}_ARTIC_medaka.json.gz")
    files_to_keep_full.append("csv_files/{sample_name}_ARTIC_medaka.csv")
    artic_results_pathlib = Path(artic_results_path)
    for filey in artic_results_pathlib.iterdir():
        # delete pangolin tree files
        if filey.is_dir():
            #protect the json for lineage analysis
            if os.path.basename(os.path.normpath(filey)) != "json_files" and os.path.basename(os.path.normpath(filey)) != "csv_files":
                rmtree(filey, ignore_errors=True)
        elif f"{filey.name}" not in files_to_keep_full:
            filey.unlink()
        elif not filey.suffix == ".gz" and filey.suffix not in [".dat", ".png"]:
            subprocess.Popen(["gzip", "-9", "-f", str(filey)]).communicate()
        elif filey.suffix == ".png":
            debug = int(get_env_variable("MT_DJANGO_DEBUG"))
            static_path = (
                f"{BASE_DIR}/artic/static/artic" if debug else f"{STATIC_ROOT}/artic"
            )
            # Copy the pngs to the artic static directory to be served
            if not Path(f"{static_path}/{artic_results_pathlib.parent.stem}").exists():
                Path(
                    f"{static_path}/{artic_results_pathlib.parent.stem}", parents=True
                ).mkdir()
            copy(
                str(filey), f"{static_path}/{artic_results_pathlib.parent.stem}",
            )

@app.task
def run_variant_command(base_results_directory, barcode_name,jm):
    """
    jm: jobmaster
    """
    re_gzip = False
    os.chdir(f"{base_results_directory}/{barcode_name}")
    MT_VoC_PATH = get_env_variable("MT_VoC_PATH")
    if (
        not Path(f"{barcode_name}.muscle.out.fasta").exists()
        and Path(f"{barcode_name}.muscle.out.fasta.gz").exists()
    ):
        logger.debug(f"Unzipping {barcode_name}.muscle.out.fasta.gz")
        subprocess.Popen(
            ["gzip", "-d", f"{barcode_name}.muscle.out.fasta.gz"]
        ).communicate()
        re_gzip = True
    cmd = [
        "bash",
        "-c",
        f"aln2type json_files csv_files {barcode_name}_ARTIC_medaka.csv MN908947.3  {barcode_name}.muscle.out.fasta {MT_VoC_PATH}/variant_definitions/variant_yaml/*.yml",
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    out, err = proc.communicate()
    if not out and err:
        logger.debug(out)
        logger.debug(err)
    else:
        logger.info("¯\_(ツ)_/¯")
        logger.info(str(out))
        logger.info(err)
    if re_gzip:
        logger.debug(f"ReGzipping {barcode_name}.muscle.out.fasta")
        subprocess.Popen(
            ["gzip", "-9", f"{barcode_name}.muscle.out.fasta"]
        ).communicate()
    try:
        VoCs_df = pd.read_csv(
            f"{barcode_name}_ARTIC_medaka.csv"
        )
        separator = ", "
        VoCs = separator.join(VoCs_df['phe-label'].to_list())

        if len(VoCs_df)> 0:
            m, created = Message.objects.get_or_create(
                recipient=jm.flowcell.owner,
                sender=jm.flowcell.owner,
                #title__startswith=f"Artic pipeline has finished for barcode {['barcode__name']}",
                flowcell=jm.flowcell,
                #defaults={
                #    "title": f"Artic pipeline has finished for barcode {barcode['barcode__name']} at {time}"
                #},
                title= f"One or more VoCs - {VoCs} - has been detected on {jm.flowcell.name} on {barcode_name}."
            )

    except FileNotFoundError as e:
        VoCs = "None Found"
        logger.error(e)
    pass

@app.task
def run_pangolin_command(base_results_directory, barcode_name):
    """

    Returns
    -------

    """
    # jm = JobMaster.objects.get(pk=job_master_pk)
    re_gzip = False
    os.chdir(f"{base_results_directory}/{barcode_name}")
    MT_CONDA_PREFIX=get_env_variable("MT_CONDA_PREFIX")
    if (
        not Path(f"{barcode_name}.consensus.fasta").exists()
        and Path(f"{barcode_name}.consensus.fasta.gz").exists()
    ):
        logger.debug(f"Unzipping {barcode_name}.consensus.fasta")
        subprocess.Popen(
            ["gzip", "-d", f"{barcode_name}.consensus.fasta.gz"]
        ).communicate()
        re_gzip = True
    #cmd = [
    #    "bash",
    #    "-c",
    #    f"source /Users/matt/.miniconda3/etc/profile.d/conda.sh && conda activate pangolin && pangolin -p --write-tree {barcode_name}.consensus.fasta",
    #]
    cmd = [
        "bash",
        "-c",
        f"source /{MT_CONDA_PREFIX} && conda activate pangolin && pangolin {barcode_name}.consensus.fasta",
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    out, err = proc.communicate()
    if not out and err:
        logger.debug(out)
        logger.debug(err)
    else:
        logger.info("¯\_(ツ)_/¯")
        logger.info(str(out))
        logger.info(err)

    if re_gzip:
        logger.debug(f"ReGzipping {barcode_name}.consensus.fasta")
        subprocess.Popen(
            ["gzip", "-9", f"{barcode_name}.consensus.fasta"]
        ).communicate()
    # if out and err:
    #    raise Exception(out)


def clear_old_data(artic_results_path, barcode_name):
    """
    TODO kind of similar to clear_unused function, could be merged?
    Clear old data of this is a rerun from the directory
    Parameters
    ----------
    artic_results_path: pathlib.PosixPath
        Path to the barcode level artic results directory
    barcode_name: str
        Name of the barcode results came from.
    Returns
    -------
    None

    """
    files_to_keep_suffix = [
        ".consensus.fasta",
        ".muscle.out.fasta",
        ".counts.dat",
        ".coverage.dat",
        ".fasta.gz",
        "lineage_report.csv.gz",
    ]
    files_to_keep = [f"{barcode_name}{filey}" for filey in files_to_keep_suffix]
    for filey in artic_results_path.iterdir():
        if filey not in files_to_keep:
            if filey.is_dir():
                rmtree(filey, ignore_errors=True)
            elif filey.is_file():
                filey.unlink()


@app.task(time_limit=1200, soft_time_limit=1200)
def run_artic_command(base_results_directory, barcode_name, job_master_pk):
    """
    Run the artic pipeline in this first pass method of dirty bash script
    Parameters
    ----------
    base_results_directory: str
        Path to the fastq file for this barcodes reads
    barcode_name: str
        The name of the barcode in string
    job_master_pk: int¯
        Primary key of the job master for this barcode
    Returns
    -------

    """
    jm = JobMaster.objects.get(pk=job_master_pk)
    jm.running = True
    jm.save()
    MT_CONDA_PREFIX=get_env_variable("MT_CONDA_PREFIX")
    # Path to barcode fastq
    fastq_path = f"{base_results_directory}/{barcode_name}/{barcode_name}.fastq"
    if not Path(fastq_path).exists():
        fastq_path += ".gz"
    # TODO get the barcode from the posix path for the sample name
    logger.info(fastq_path)
    scheme_dir = get_env_variable("MT_ARTIC_SCEHEME_DIR")
    artic_env = get_env_variable("MT_ARTIC_ENV")
    os.chdir(f"{base_results_directory}/{barcode_name}")
    cmd = [
        "bash",
        "-c",
        f"source {MT_CONDA_PREFIX} && conda activate {artic_env} && artic minion --medaka --normalise 200 --threads 4 --scheme-directory {scheme_dir} --read-file {fastq_path} nCoV-2019/V3 {barcode_name}",
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    if not out and err:
        logger.debug("We are done here. Artic out!")
        logger.debug(out)
        logger.debug(err)
    else:
        logger.debug("¯\_(ツ)_/¯")
        logger.warning(str(out))
        logger.warning(err)

    # Now run the aln2type artic screen.

    run_variant_command(base_results_directory, barcode_name, jm)

    run_pangolin_command(base_results_directory, barcode_name)
    # Update the Barcode Metadata to show the task has been run on this barcode

    ArticBarcodeMetadata.objects.filter(
        flowcell=jm.flowcell, job_master__job_type_id=16, barcode__name=barcode_name
    ).update(has_finished=True, marked_for_rerun=False)

    clear_unused_artic_files(
        f"{base_results_directory}/{barcode_name}", barcode_name, jm.flowcell.id
    )
    jm = JobMaster.objects.get(pk=job_master_pk)
    jm.running = False
    jm.complete = True
    jm.save()


def save_artic_command_job_masters(flowcell, barcode_name, reference_info, run_id):
    """
    Save the JobMasters that we need to manage the queue of tasks.
    Parameters
    ----------
    flowcell: minknow_data.models.Flowcell
        Flowcell these runs are on
    barcode_name: str
        Barcode name for these reads
    reference_info: reference.models.ReferenceInfo
        Covid Reference
    run_id: int
        The primary key of the run
    Returns
    -------

    """
    # Don't automatically run unclassified barcode
    if barcode_name == "unclassified":
        return
    job_type = JobType.objects.get(name="Run Artic")
    # TODO potential bug here where the barcode name on the reads is not the barcode name we save and alos multiple runs. if force unique is false
    barcode_object = Barcode.objects.get(run_id=run_id, name=barcode_name)
    job_master, created = JobMaster.objects.get_or_create(
        job_type=job_type,
        reference=reference_info,
        barcode=barcode_object,
        flowcell=flowcell,
    )
    if not created:
        logger.info("Firing artic due to condition being met again")
        # todo the below could be a problem
        job_master.running = False
        job_master.complete = False
        job_master.save()
    return


def make_barcoded_directories(barcodes, results_dir):
    """

    Parameters
    ----------
    barcodes: list of str
        List of barcode names
    results_dir: pathlib.PosixPath
        Path to the base results directory

    Returns
    -------
    None

    """
    for barcode in barcodes:
        path = results_dir / barcode
        if not path.exists():
            Path.mkdir(path, parents=True)
            logger.debug(f"Created directory {path} for barcode {barcode}.")
        else:
            logger.debug(f"Directory {path} already exists.")


def add_barcode_to_tofire_file(barcode_name, directory):
    """
    Add a barcode to our file of barcodes we have fired the file for.
    Parameters
    ----------
    barcode_name: str
        Barcode name that we have enough coverage for to fire pipeline.
    directory: pathlib.PosixPath
        Absolute path to the results directory that we are storing this file in.
    Returns
    -------
    exit_code: int
        If 0 successfully added, else returns 1 and the error.
    """
    path = directory / "barcodes_to_fire.txt"
    with open(path, "a") as fh:
        try:
            fh.write(f"{barcode_name}\n")
            return 0
        except Exception as e:
            logger.error(str(e))
            return 1, e


def fetch_barcode_to_fire_list(directory):
    """

    Parameters
    ----------
    directory: pathlib.PosixPath
        the absolute path to the top level of th results directory of this task
    Returns
    -------
    barcode: list of str
        List of barcode to fire.

    """

    path = directory / "barcodes_to_fire.txt"
    try:
        with open(path, "r") as fh:
            barcodes = fh.read().split("\n")
    except FileNotFoundError as e:
        logger.debug("File not found. First iteration? Continuing...")
        return []

    return barcodes


def populate_reference_count(
    reference_count_dict, fp, master_reference_count_path,
):
    """
    Populate the reference count dict so we can write it out with counts
    Parameters
    ----------
    reference_count_dict: dict
        Empty dictionary
    fp: file
        Filepath to the reference file
    master_reference_count_path: pathlib.PosixPath
        Path to the reference fasta on the server

    Returns
    -------
    reference_count_dict: dict
        The time of your life

    """
    # # Use the reference to populate the priors dict
    # If priors_dict and reference_count_dict for this dataset ALREADY EXISTS load them from archive

    if master_reference_count_path.exists():
        with open(master_reference_count_path, "rb") as f2:
            reference_count_dict = pickle.load(f2)
    else:
        for desc, name, seq, qual in readfq(fp):
            # Call multi array results to create a structured array of zeros the length of the reference,
            #  with 9 fields, A,C,G,T,D,I,IC,M,U
            reference_count_dict[name] = multi_array_results(len(seq))
            # create a dictionary under the priors dict, keyed to the reference name

        # Save reference_count_dict for next use
        with open(master_reference_count_path, "wb") as fh2:
            pickle.dump(reference_count_dict, fh2, pickle.HIGHEST_PROTOCOL)

    return reference_count_dict


def replicate_counts_array_barcoded(barcode_names, counts_dict):
    """
    Replicate the empty counts numpy array for each barcode
    Parameters
    ----------
    barcode_names: list of str
        the name of the barcode
    counts_dict: dict
        Dict keyed to chromosome/contig name. Value of the array for the counts.
        Array length of the number of bases in reference.

    Returns
    -------
    barcoded_counts_dict: dict
        Dictionary to store results. Keyed barcode name to a value of the reference count dictionary.
    None
    """
    barcoded_counts_dict = {barcode: deepcopy(counts_dict) for barcode in barcode_names}
    return barcoded_counts_dict


def save_artic_barcode_metadata_info(
    coverage,
    flowcell,
    job_master,
    run_id,
    barcode_name,
    has_sufficient_coverage=False,
    projected_to_finish=False,
):
    """
    Save the artic metadata info for this barcode.
    Parameters
    ----------
    coverage: np.ndarray
        Array of coverage, element for each base.

    flowcell: minknow_data.models.Flowcell
        Flowcell Django ORM object for the flowcell that produced these reads.

    job_master: reads.models.JobMaster
        JobMaster Django ORM object for the JobMaster representing this running task.

    run_id: int
        Primary key of the run record in the database.

    barcode_name: str
        The name of the barcode.

    has_sufficient_coverage: bool
        If True, the barcode has sufficient coverage to fire the artic command

    projected_to_finish: bool
        If the barcode is projected to finish by our algorithm

    Returns
    -------
    None
    """
    try:
        barcode_numbers = FlowcellSummaryBarcode.objects.filter(
            flowcell=flowcell, barcode_name=barcode_name, status="True"
        )
        all_numbers = FlowcellSummaryBarcode.objects.filter(
            flowcell=flowcell, barcode_name="All reads"
        ).values_list("read_count", flat=True)
    except FlowcellSummaryBarcode.DoesNotExist as e:
        logger.warning(e)
        return
    barcode_numbers_total = 0
    for barcode_number in barcode_numbers:
        barcode_numbers_total += barcode_number.read_count
    total_reads = 0
    for all_number in all_numbers:
        total_reads += all_number
    # Except zero division in case no flowcell summary barcode entries yet
    try:
        proportion = barcode_numbers_total / total_reads * 100
    except ZeroDivisionError:
        proportion = 0
    # #### Save coverage data for the barcode
    try:
        barcode = Barcode.objects.get(run_id=run_id, name=barcode_name)
    except Barcode.DoesNotExist as e:
        logger.warning(
            f"Barcode {barcode_name} not found for run_id {run_id}. Skipping..."
        )
        return
    barcodes_minimum_coverage = round(coverage.min(), 2)
    barcodes_maximum_coverage = round(coverage.max(), 2)
    barcodes_variance_coverage = round(coverage.var(), 2)
    barcodes_average_coverage = round(coverage.mean(), 2)
    defaults = {
            "average_coverage": barcodes_average_coverage,
            "maximum_coverage": barcodes_maximum_coverage,
            "minimum_coverage": barcodes_minimum_coverage,
            "variance_coverage": barcodes_variance_coverage,
            "percentage_of_reads_in_barcode": round(proportion, 2),
            "projected_to_finish": projected_to_finish
        }
    if has_sufficient_coverage:
        defaults.update({"has_sufficient_coverage": has_sufficient_coverage})
    orm_object, created = ArticBarcodeMetadata.objects.update_or_create(
        flowcell=flowcell,
        job_master=job_master,
        barcode=barcode,
        defaults=defaults
    )


@app.task(base=MyTask)
def run_artic_pipeline(task_id, streamed_reads=None):
    """
    Run the artic pipeline on a flowcells amount of reads
    Parameters
    ----------
    task_id: int
        The pk of the task record in the database
    streamed_reads: list of dict
        List of reads in dictionary form. Default None.

    Returns
    -------

    """
    # The location of the mimimap2 executable
    minimap2 = getattr(settings, "MINIMAP2", None)
    task = JobMaster.objects.get(pk=task_id)
    if not task.reference:
        raise ValueError("Missing Reference file. Please sort out.")
    flowcell = task.flowcell
    avg_read_length = flowcell.average_read_length
    if avg_read_length == 0:
        logger.error(
            f"Average read length is zero Defaulting to 450, but this is an error."
        )
        avg_read_length = 450
    if not streamed_reads:
        # The chunk size of reads from this flowcell to fetch from the database
        # aim for 50 megabases
        desired_yield = 50 * 1000000
        chunk_size = round(desired_yield / avg_read_length)
        logger.debug(f"Fetching reads in chunks of {chunk_size} for alignment.")
        fasta_objects = FastqRead.objects.filter(
            flowcell_id=flowcell.id, id__gt=task.last_read
        )[:chunk_size]
        fasta_df_barcode = pd.DataFrame().from_records(
            fasta_objects.values("read_id", "barcode_name", "sequence", "id", "run_id")
        )
        if fasta_objects:
            last_read = fasta_df_barcode.tail(1).iloc[0].id
        if fasta_df_barcode.shape[0] == 0:
            logger.debug("No fastq found. Skipping this iteration.")
            task.running = False
            task.save()
            return
    else:
        last_read = task.last_read
        fasta_df_barcode = pd.DataFrame(streamed_reads)
        fasta_df_barcode = fasta_df_barcode[fasta_df_barcode["is_pass"]]
        fasta_objects = streamed_reads
    read_count = fasta_df_barcode.shape[0]
    logger.debug(f"Fetched reads.")
    logger.debug(fasta_df_barcode.shape)
    min_read_length = 400
    max_read_length = 700
    fasta_df_barcode["sequence_length"] = fasta_df_barcode["sequence"].str.len()
    fasta_df_barcode = fasta_df_barcode[
        fasta_df_barcode["sequence_length"].between(
            min_read_length, max_read_length, inclusive=True
        )
    ]
    logger.debug(
        f"Reads after filtering for read lengths between"
        f" {min_read_length} - {max_read_length}: {fasta_df_barcode.shape[0]}"
    )
    # ##Check if read_count is > 0 otherwise exit
    ## #We might want to check the minimum number of reads or the time
    # since the last read was entered to do a smaller chunk?
    if read_count > 0:
        if not minimap2:
            logger.error("Can not find minimap2 executable - stopping task.")
            # This won't reset the task?
            return
        reference_info = task.reference
        # Chromosome dict contains all the chromsome names keyed to all the info we have about that chromosome
        chromdict = dict()
        # Get the chromosome we have from this reference
        chromosomes = reference_info.reference_lines.all()
        # For each chromosome
        for chromosome in chromosomes:
            chromdict[chromosome.line_name] = chromosome

        # Create the series that contains the read_id and sequence as a correctly formatted fasta string
        fasta_df_barcode["fasta"] = (
            ">" + fasta_df_barcode["read_id"] + "\n" + fasta_df_barcode["sequence"]
        )
        barcodes = fasta_df_barcode["barcode_name"].unique()
        # Create one string formatted correctly for fasta with input data
        fasta_data = "\n".join(list(fasta_df_barcode["fasta"]))
        if streamed_reads:
            fastq_dict = {fasta.get("read_id"): fasta for fasta in fasta_objects}
        else:
            # dict where we will be storing the fastq objects keyed to read ID
            fastq_dict = {fasta.read_id: fasta for fasta in fasta_objects}
        # create the command we are calling minimap with
        cmd = "{} -x map-ont -t 1 --secondary=no -c --MD {} -".format(
            minimap2, reference_info.file_location.path.strip()
        )
        logger.debug(
            "Flowcell id: {} - Calling minimap Artic - {}".format(flowcell.id, cmd)
        )
        # Setup minimap call
        proc = subprocess.Popen(
            cmd.split(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdin=subprocess.PIPE,
            universal_newlines=True,
        )
        # SubProcess communicate actually calls the function
        (out, err) = proc.communicate(input=fasta_data)
        if err:
            logger.error(err)
        logger.debug(
            "Flowcell id: {} - Finished minimap artic - {}".format(flowcell.id, cmd)
        )
        # decode output from byte string
        paf = out
        if paf:
            reference_count_dict = {}
            chromosomes_seen_now = set()
            # Dictionary to store reads that have mapped under the correct barcode_name
            barcode_sorted_fastq_cache = defaultdict(list)
            reference_path = reference_info.file_location.path
            base_result_dir_path = make_results_directory_artic(flowcell.id, task.id)
            # We can save time by writing the empty numpy data struct and reloading it
            master_reference_count_path = Path(
                f"{base_result_dir_path}/{reference_info.name}_counts.pickle"
            )
            if reference_info.file_name.endswith(".gz"):
                with gzip.open(reference_path, "rt") as fp:
                    reference_count_dict = populate_reference_count(
                        reference_count_dict, fp, master_reference_count_path
                    )
            else:
                with open(reference_path, "r") as fp:
                    reference_count_dict = populate_reference_count(
                        reference_count_dict, fp, master_reference_count_path
                    )
            run_id = fasta_df_barcode["run_id"].unique()[0]
            barcodes = np.unique(fasta_df_barcode["barcode_name"])
            for barcode_name in barcodes:
                barcode_orm = Barcode.objects.get(run_id=run_id, name=barcode_name)
                orm_object, created = ArticBarcodeMetadata.objects.get_or_create(
                    flowcell=flowcell, job_master=task, barcode=barcode_orm
                )
            logger.debug(f"Barcodes in this iteration are {barcodes}")
            barcoded_counts_dict = replicate_counts_array_barcoded(
                barcodes, reference_count_dict
            )
            paf_summary_cov_dict = {}
            make_barcoded_directories(barcodes, base_result_dir_path)
            logger.debug(f"Parsing paf file. Please wait.")
            barcodes_with_mapped_results = set()
            for i, record in enumerate(
                parse_PAF(
                    StringIO(paf),
                    fields=[
                        "qsn",
                        "qsl",
                        "qs",
                        "qe",
                        "rs",
                        "tsn",
                        "tsl",
                        "ts",
                        "te",
                        "nrm",
                        "abl",
                        "mq",
                        "tags",
                    ],
                )
            ):
                fastq_read = fastq_dict[record.qsn]
                if streamed_reads:
                    barcode_of_read = fastq_read.get("barcode_name", "Failed")
                    read_id = fastq_read.get("read_id")
                    sequence = fastq_read.get("sequence")
                    read_type_pk = fastq_read.get("type", 1)
                else:
                    barcode_of_read = fastq_read.barcode_name
                    read_id = fastq_read.read_id
                    sequence = fastq_read.sequence
                    read_type_pk = fastq_read.type.id
                barcodes_with_mapped_results.add(barcode_of_read)
                barcode_sorted_fastq_cache[barcode_of_read].append(
                    f">{read_id}\n{sequence}"
                )
                chromosome = chromdict[record.tsn]
                #########################################################
                # ########## Now for the parse cigar string ########### #
                #########################################################
                # Results array
                mapping_start = record.ts
                # The base at which the mapping ends
                mapping_end = record.te
                # The length of the mapping
                mapping_length = mapping_end - mapping_start
                if barcode_of_read not in paf_summary_cov_dict:
                    paf_summary_cov_dict[barcode_of_read] = {
                        "job_master": task,
                        "barcode_name": barcode_of_read,
                        "chromosome": chromosome,
                        "reference_line_length": chromosome.chromosome_length,
                        "yield": 0,
                        "read_count": 0,
                        "read_type_id": read_type_pk,
                    }
                paf_summary_cov = paf_summary_cov_dict[barcode_of_read]
                paf_summary_cov["yield"] += mapping_length
                paf_summary_cov["read_count"] += 1
                chromosomes_seen_now.add(chromosome.line_name)
                # Parse the cigar string, on the path file record, returns a dictionary of bases, with numpy arrays
                d = parse_cigar(record.tags.get("cg", None), sequence, mapping_length,)
                # Adds the mismatch count to the temporary dict, parsing the MD string
                d["M"] = np.fromiter(parse_MD(record.tags.get("MD", None)), np.uint16)
                # Loop through the dictionary of results, through each base
                for base in d:
                    # add the counts for that base to the array that exists in the reference count dict for that base,
                    # at the correct location
                    barcoded_counts_dict[barcode_of_read][record.tsn][base][
                        mapping_start:mapping_end
                    ] += d[base].astype(reference_count_dict[record.tsn][base].dtype)
                # TODO this needs to be reset, here we set the U to True to indicate that a count at this position
                # ToDO has been changed this iteration, might offer speed up if we do this?
                # reference_count_dict[record.tsn]["U"][
                # mapping_start:mapping_end
                # ] = True
            # iterate the paf_summary_cov objects
            read_counts_dict = {}
            for paf_summary_barcode, paf_summary_cov in paf_summary_cov_dict.items():
                logger.debug(f"Creating PafSummaryCov objects for {paf_summary_barcode}")
                paf_summary_cov_orm, created = PafSummaryCov.objects.get_or_create(
                    job_master=task,
                    barcode_name=paf_summary_cov["barcode_name"],
                    chromosome=paf_summary_cov["chromosome"],
                    reference_line_length=paf_summary_cov["reference_line_length"],
                    read_type_id=paf_summary_cov["read_type_id"],
                )
                paf_summary_cov_orm.total_yield += paf_summary_cov["yield"]
                paf_summary_cov_orm.read_count += paf_summary_cov["read_count"]
                read_counts_dict[
                    paf_summary_cov["barcode_name"]
                ] = paf_summary_cov_orm.read_count
                paf_summary_cov_orm.coverage = round(
                    paf_summary_cov_orm.total_yield
                    / paf_summary_cov_orm.reference_line_length,
                    2,
                )
                paf_summary_cov_orm.average_read_length = round(
                    paf_summary_cov_orm.total_yield / paf_summary_cov_orm.read_count
                )
                paf_summary_cov_orm.save()
            # TODO only ever see one chromosome, so we can remove for loop?
            barcodes_already_fired = fetch_barcode_to_fire_list(base_result_dir_path)
            num_amplicons, amplicon_band_coords = get_amplicon_infos()
            time_stamp = datetime.now()
            df_old = get_or_create_metrics_df(base_result_dir_path)
            df_new_dict = {}
            # TODO split into own function
            for chrom_key in chromosomes_seen_now:
                # Write out the coverage
                for barcode_name in barcodes_with_mapped_results:
                    barcode_sorted_fastq_path = Path(
                        f"{base_result_dir_path}/{barcode_name}/{barcode_name}.fastq"
                    )
                    # TODO NOTE REMOVED CHROMOSOME FROM PATH SO WILL NOT WORK FOR MULTI CHROMOSOME DATA
                    coverage_path = Path(
                        f"{base_result_dir_path}/{barcode_name}/{barcode_name}.coverage.dat"
                    )
                    counts_path = Path(
                        f"{base_result_dir_path}/{barcode_name}/{barcode_name}.counts.dat"
                    )
                    ############################################################################
                    # ###################### read in the Counts ############################ #
                    # ##################### Combine extant values ############################ #
                    ############################################################################
                    max_values = 0
                    if counts_path.exists():
                        # Read the old coverage counts into memory
                        with open(counts_path, "rb") as fh:
                            old_counts_array = np.fromfile(
                                fh,
                                dtype=[
                                    ("A", np.uint16),
                                    ("C", np.uint16),
                                    ("G", np.uint16),
                                    ("T", np.uint16),
                                    ("D", np.uint16),
                                    ("I", np.uint16),
                                    ("IC", np.uint16),
                                    ("M", np.uint16),
                                    ("U", np.bool),
                                ],
                            )
                        # Add the sum of each position counts to the old array, by each base type
                        for name in barcoded_counts_dict[barcode_name][
                            chrom_key
                        ].dtype.names:
                            # Don't update U as True + True = 2, which would ruin everything
                            if name != "U":
                                barcoded_counts_dict[barcode_name][chrom_key][
                                    name
                                ] += old_counts_array[name]
                                a = barcoded_counts_dict[barcode_name][chrom_key][name].max()
                                if a > max_values:
                                    max_values = a
                    coverage = np.sum(
                        np.array(
                            barcoded_counts_dict[barcode_name][chrom_key][
                                ["A", "C", "G", "T"]
                            ].tolist(),
                            dtype=np.int16,
                        ),
                        axis=1,
                        dtype=np.uint16,
                    )
                    # ################################################################ #
                    # ################# Write out coverage array ################ #
                    ####################################################################
                    # Add the max value for each of the arrays
                    with open(coverage_path, "wb") as fh:
                        fh.write(coverage.data)
                    # ################################################################ #
                    # ################### Write out counts arrays ######################
                    ####################################################################
                    with open(counts_path, "wb") as fhc:
                        fhc.write(barcoded_counts_dict[barcode_name][chrom_key].data)
                    with open(barcode_sorted_fastq_path, "a") as fh:
                        fh.write("\n".join(barcode_sorted_fastq_cache[barcode_name]))
                    logger.debug(
                        f"Avg Coverage at iteration {task.iteration_count} for barcode_name {barcode_name} is {coverage.mean()}"
                    )
                    # get artic fire conditions
                    afcs = ArticFireConditions.objects.filter(flowcell=flowcell)
                    has_sufficient_coverage = False
                    amplicon_stats = get_amplicon_stats(
                        amplicon_band_coords=amplicon_band_coords,
                        num_amplicons=num_amplicons,
                        flowcell_id=flowcell.id,
                        barcode_name=barcode_name,
                        time_stamp=time_stamp,
                        read_count=read_counts_dict[barcode_name],
                        coverage_array=coverage,
                    )
                    if any(
                        [
                            check_afc_values_met(
                                afc,
                                amplicon_stats.amplicon_coverage_medians,
                                num_amplicons,
                                run_id,
                                barcode_name,
                                task
                            )
                            for afc in afcs
                        ]
                    ):
                        if barcode_name not in barcodes_already_fired:
                            add_barcode_to_tofire_file(barcode_name, base_result_dir_path)
                        save_artic_command_job_masters(
                            flowcell, barcode_name, reference_info, run_id
                        )
                        has_sufficient_coverage = True
                    # save the artic per barcode_name metadata
                    projected_to_finish = predict_barcode_will_finish(
                        amplicon_stats.amplicon_coverage_medians,
                        flowcell.number_barcodes,
                        task.read_count + read_count,
                    )
                    save_artic_barcode_metadata_info(
                        coverage,
                        flowcell,
                        task,
                        run_id,
                        barcode_name,
                        has_sufficient_coverage,
                        projected_to_finish
                    )
                    amp_stat_dict = amplicon_stats._asdict()
                    amp_stat_dict["projected_to_finish"] = projected_to_finish
                    amp_stat_dict["total_read_count"] = task.read_count + read_count
                    df_new_dict.update(
                        {
                            barcode_name: amp_stat_dict
                        }
                    )
            df_new = pd.DataFrame.from_dict(df_new_dict, orient="index")
            write_out_artic_metrics(
                df_old=df_old, df_new=df_new, dir_to_write_to=base_result_dir_path
            )
    task.last_read = last_read
    task.iteration_count += 1
    logger.debug("Finishing this batch of reads.")
    task.read_count += read_count
    task.running = False
    task.save()
