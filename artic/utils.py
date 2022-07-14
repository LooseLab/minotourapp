"""
Useful functions that don't belong in views or tasks
"""
import json
import math
import os
import subprocess
import tarfile
from collections import namedtuple
from datetime import datetime
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd
import toyplot.pdf
import toyplot.svg
import toytree as toytree
from celery.utils.log import get_task_logger
from django.http import HttpResponse
from django.shortcuts import get_object_or_404, render
from weasyprint import CSS, HTML

from alignment.models import PafSummaryCov
from artic.models import ArticBarcodeMetadata
from minknow_data.models import Flowcell
from minotourapp.celery import app
from minotourapp.settings import BASE_DIR
from minotourapp.utils import get_env_variable
from reads.models import JobMaster

# Colour palette of amplicon bands
colour_palette = ["#ffd9dc", "#ffefdc", "#ffffbc", "#dcffe4", "#bae1ff"]
logger = get_task_logger(__name__)


def get_sequencing_stats_pdf():
    """
    Get the stats for a sequencing report into a pdf
    """
    pass


def unique_amplicon_coordinates(scheme_bed_file):
    """
    Filter out overlaps between amplicon scheme files
    Parameters
    ----------
    scheme_bed_file: Path
        Full path to the scheme file
    Returns
    -------
    numpy.ndarray

    """
    df = pd.read_csv(
        scheme_bed_file,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["chromosome", "start", "end", "name"],
    )
    df["primer_position"] = df["name"].str.split("_").str[1]
    df = df.set_index("primer_position")
    df[["primer_start", "primer_end"]] = df.groupby("primer_position").agg(
        {"start": np.min, "end": np.max}
    )
    df = df.loc[~df.index.duplicated(keep="last")]
    df["primer_end"] = (
        (df["primer_start"].shift(-1) - 1).fillna(df["primer_end"])
    ).astype(int)
    df["pair_end"] = (
        (df["primer_start"].shift(-1) - 1).fillna(df["primer_end"])
    ).astype(int)
    df["previous_end"] = df["end"].shift()
    unqiue_amplicon_coords = np.column_stack(
        (
            np.where(
                df["primer_start"] < df["previous_end"],
                df["previous_end"],
                df["primer_start"],
            ),
            df["primer_end"].values,
        )
    )
    unqiue_amplicon_coords = unqiue_amplicon_coords.astype(int)
    return unqiue_amplicon_coords


def get_artic_run_stats(pk, svg_data, request, task, logger):
    """
    Generate a one page pdf on the artic run that we generated
    Parameters
    ----------
    pk: int
        The flowcell primary key
    svg_data: dict
        Dictionary containing sgv strings for each chart to be drawn onto PDF
    request: rest_framework.request.Request
        AJAX request instance
    task: reads.models.JobMaster
        The job master ORM instance for this artic task
    logger: logging.logger
        Logging instance for celery task
    Returns
    -------
    BytesIO
        Bytes of the pdf we have generated
    """
    flowcell, artic_results_path, jm_id, _ = quick_get_artic_results_directory(pk)
    artic_task = JobMaster.objects.get(pk=jm_id)
    queryset = PafSummaryCov.objects.filter(job_master=artic_task).values(
        "barcode_name",
        "chromosome__line_name",
        "average_read_length",
        "coverage",
        "job_master__id",
        "read_count",
        "total_yield",
        "reference_line_length",
    )
    artic_metadata = {
        d.barcode.name: d
        for d in ArticBarcodeMetadata.objects.filter(
            flowcell_id=pk, job_master=artic_task
        )
    }
    # dictionaries change in place
    scheme = task.primer_scheme.scheme_species
    scheme_version = task.primer_scheme.scheme_version
    scheme_dir = task.primer_scheme.scheme_directory
    amplicon_band_coords, colours = get_amplicon_band_data(
        scheme, scheme_version, scheme_dir
    )
    num_amplicons = len(amplicon_band_coords)
    time_stamp = datetime.now()
    for paf_summary_cov in queryset:
        barcode_name = paf_summary_cov["barcode_name"]
        amplicon_stats = get_amplicon_stats(
            amplicon_band_coords,
            num_amplicons,
            pk,
            barcode_name,
            time_stamp,
            paf_summary_cov.get("read_count", "Unknown"),
        )
        paf_summary_cov["has_finished"] = artic_metadata[barcode_name].has_finished

        # get the lineage if it's finished
        if paf_summary_cov["has_finished"]:
            try:
                lineage = pd.read_csv(
                    artic_results_path / barcode_name / "lineage_report.csv.gz"
                )["lineage"][0]
            except FileNotFoundError:
                lineage = "Unknown"

            try:
                VoCs_df = pd.read_csv(
                    artic_results_path
                    / barcode_name
                    / f"{barcode_name}_ARTIC_medaka.csv.gz"
                )
                VoCs = VoCs_df["phe-label"]
            except FileNotFoundError:
                VoCs = "None Found"
        else:
            lineage = "Currently unknown"
            VoCs = "Not Tested"
        if barcode_name == "unclassified":
            paf_summary_cov["has_sufficient_coverage"] = "ignore"
        else:
            paf_summary_cov["has_sufficient_coverage"] = artic_metadata[
                barcode_name
            ].has_sufficient_coverage
        paf_summary_cov["VoC_Warn"] = VoCs
        amp_stats_dict = amplicon_stats._asdict()
        amp_stats_dict.pop("amplicon_coverage_medians")
        amp_stats_dict.pop("amplicon_coverage_means")
        paf_summary_cov.update(amp_stats_dict)
        paf_summary_cov["lineage"] = lineage
        paf_summary_cov["projected_to_finish"] = artic_metadata[
            barcode_name
        ].projected_to_finish
    svg_data["sample_name"] = flowcell.name
    l = sorted(list(queryset), key=lambda x: x["barcode_name"])
    svg_path = artic_results_path / "snp_plot.svg"
    tree_path = artic_results_path / "iqtree_.treefile"
    logger.info("treeeee")
    svg_data["treesnstuff"] = False
    if tree_path.exists():
        with open(tree_path, "r",) as fh:
            trey = toytree.tree(fh.read(), tree_format=0)
            canvas, axes, mark = trey.draw(
                width=1000,
                height=700,
                node_sizes=12,
                node_style={"fill": "green", "stroke": "black", "stroke-width": 0.75,},
            )
            tree_path_png = f"/tmp/{flowcell.id}_tree-plot.pdf"
            toyplot.pdf.render(canvas, tree_path_png)
        # svg_data["treesnstuff"] = True
        # svg_data["tree_path_svg"] = f"file:///{tree_path_png}"
    logger.info("no treeeee")

    svg_data["overall_results"] = l
    # print(svg_data)
    artic_full_report_bytes = BytesIO(
        HTML(
            string=render(
                request, "artic-report.html", context={"data": svg_data}
            ).getvalue()
        ).write_pdf(
            None,
            stylesheets=[
                CSS(f"{BASE_DIR}/web/static/web/css/artic-report.css"),
                CSS(f"{BASE_DIR}/web/static/web/libraries/bootstrap-4.5.0-dist/css/bootstrap.css"),
            ],
        )
    )
    return artic_full_report_bytes


def get_all_results(artic_results_dir, flowcell, selected_barcode, chosen):
    """
    
    Parameters
    ----------
    artic_results_dir: pathlib.PosixPath
        The path to the artic results directory
    flowcell: minknow_data.models.Flowcell
        The ORM object of the flowcell
    selected_barcode: str
        The name of the selected barcode, empty string if downloading all barcodes
    chosen: list
        List of files that were chosen for download
    Returns
    -------
    list
        List of completed barcodes for this flowcell
    """
    if not selected_barcode:
        queryset = ArticBarcodeMetadata.objects.filter(
            flowcell=flowcell, has_finished=True
        )
        finished_barcodes = [q.barcode.name for q in queryset]
    else:
        finished_barcodes = [selected_barcode]
    results_files = {
        "consensus": [
            f"{barcode_name}/{barcode_name}.consensus.fasta"
            for barcode_name in finished_barcodes
        ],
        "fail-vcf": [
            f"{barcode_name}/{barcode_name}.fail.vcf.gz"
            for barcode_name in finished_barcodes
        ],
        "pass-vcf": [
            f"{barcode_name}/{barcode_name}.pass.vcf.gz"
            for barcode_name in finished_barcodes
        ],
        "input-fasta": [
            f"{barcode_name}/{barcode_name}.fastq"
            for barcode_name in finished_barcodes
        ],
        "pangolin-lineages": [
            f"{barcode_name}/lineage_report.csv"
            for barcode_name in finished_barcodes
        ],
        "sorted-bam": [
            f"{barcode_name}/{barcode_name}.sorted.bam"
            for barcode_name in finished_barcodes
        ],
        "sorted-bam-bai": [
            f"{barcode_name}/{barcode_name}.sorted.bam.bai"
            for barcode_name in finished_barcodes
        ],
    }
    chosen_files = [results_files[key] for key in chosen]

    # change into the directory
    os.chdir(artic_results_dir)
    results_file = artic_results_dir / f"results_artic_{flowcell.name}.tar.gz"
    with tarfile.open(results_file, "w:gz") as tar:
        for filey in chosen_files:
            try:
                for barcode_file in filey:
                    tar.add(barcode_file)
            except FileNotFoundError as e:
                print(f"file not found {repr(e)}")
                try:
                    for barcode_file in filey:
                        tar.add(f"{barcode_file}.gz")
                except FileNotFoundError as e:
                    print(f"file not found {repr(e)}")

    with open(results_file, "rb") as fh:
        response = HttpResponse(fh.read(), content_type="application/gzip")
        response[
            "Content-Disposition"
        ] = f"attachment; filename=results_artic_{flowcell.name}.tar.gz"
        response["x-file-name"] = f"results_artic_{flowcell.name}.tar.gz"
        return response


def get_amplicon_band_data(scheme, scheme_version, scheme_dir, nooverlap=True):
    """
    Retrieve coordinates on reference for amplicon bands, and a colour scheme for any amplicon pools
    Parameters
    ----------
    scheme: str
        Artic primer scheme to get data from
    scheme_version: str
        The version, if any
    scheme_dir: str
        The scheme directory path
    Returns
    -------
    (list of string, dict)
        A list of coordinates, start and stop on x axis, and a colour scheme lookup dictionary for amplicon pools.

    """
    json_file_path = get_amplicon_json_file_path(
        scheme, scheme_version, Path(scheme_dir)
    )
    with open(json_file_path, "r") as fh:
        amplicon_bands = json.load(fh)
    # Get data
    if nooverlap:
        amplicon_band_coords = json.loads(amplicon_bands["unique_amplicons"])
    else:
        amplicon_band_coords = json.loads(amplicon_bands["amplicons"])
    colours = {
        amplicon_bands["pools"][index]: colour_palette[index]
        for index in range(len(amplicon_bands["pools"]))
    }
    return amplicon_band_coords, colours


def get_amplicon_json_file_path(scheme, scheme_version, scheme_dir):
    """
    Get the file path to the amplicon primer JSON file
    Parameters
    ----------
    scheme: Path
        Artic primer scheme to get data from
    scheme_version: str
        The version, if any
    scheme_dir: str
        The scheme directory
    Returns
    -------
    pathlib.PosixPath
        File path to the JSON primer file.

    """
    artic_results_dir_primers = (
        Path(get_env_variable("MT_ARTIC_RESULTS_DIR")) / "artic/" / "primers_files"
    )
    if not (Path(get_env_variable("MT_ARTIC_RESULTS_DIR")) / "artic/").exists():
        (Path(get_env_variable("MT_ARTIC_RESULTS_DIR")) / "artic/").mkdir()
    # Needs to be dynamic here
    bed_file = f"{scheme}.scheme.bed"
    json_file = f"{scheme}_{scheme_version}.primers.json"
    full_path_to_bed_file = scheme_dir / scheme / scheme_version / bed_file
    if not (artic_results_dir_primers / json_file).exists():
        json_file_path = convert_amplicon_bed_file_to_json(
            full_path_to_bed_file, json_file, artic_results_dir_primers
        )
    else:
        json_file_path = artic_results_dir_primers / json_file
    return json_file_path


def convert_amplicon_bed_file_to_json(filepath, json_file, artic_results_primer_dir):
    """
    TODO nowhere near dynamic enough
    Convert an amplicon primer scheme bed file to JSON listing the primer start and end
    Parameters
    ----------
    filepath: Path
        Absolute file path to the primer scheme bed file
    json_file: str
        Name to give the JSON file created from the bed file
    artic_results_primer_dir: Path
        Absolute path to the directory to store amplicon JSON in

    Returns
    -------
    pathlib.PosixPath
        File path to JSON file, stored in Artic Results Dir
    """
    # TODO not dynamic, how make dynamic, very hardcoded
    # TODO Dropdown on task start for scheme dir?
    df = pd.read_csv(filepath, sep="\t", header=None, usecols=[0,1,2,3])
    df = df[df.columns[~df.isnull().all()]]
    df["primer_number"] = pd.to_numeric(df[3].str.split("_").str[1])
    df = df.set_index("primer_number")
    df[["primer_start", "primer_end"]] = df.groupby("primer_number").agg(
        {1: np.min, 2: np.max}
    )
    df = df.reset_index()
    df = df.set_index(["primer_start", "primer_end"])
    df = df.loc[~df.index.duplicated(keep="first")]
    df = df.reset_index()
    df[["unique_primer_start", "unique_primer_end"]] = unique_amplicon_coordinates(filepath)
    json_data = {
        "amplicons": df[["primer_start", "primer_end", 4]].to_json(orient="values"),
        "unique_amplicons":df[["unique_primer_start", "unique_primer_end", 4]].to_json(orient="values"),
        "name": f"{df[0].unique()[0]}_primer_scheme",
        "pools": df[4].unique().tolist(),
    }
    # Check bed file results dir exists
    if not artic_results_primer_dir.exists():
        artic_results_primer_dir.mkdir()
    with open(f"{artic_results_primer_dir}/{json_file}", "w") as fh:
        json.dump(json_data, fh)
    return Path(f"{artic_results_primer_dir}/{json_file}")


def quick_get_artic_results_directory(flowcell_id, barcodeName="", check=False):
    """

    Parameters
    ----------
    flowcell_id : int
        Primary key of flowcell record that we are looking
    barcodeName : str optional
        The barcode that the data is from. If provided coverage path is also returned.
    check : bool
        If True, check whether there ia a results directory.

    Returns
    -------

    (flowcell: reads.models.Flowcell, artic_results_path: pathlib.PosixPath, artic_task.id : int)



    """
    flowcell = Flowcell.objects.get(pk=flowcell_id)
    artic_task = get_object_or_404(JobMaster, flowcell=flowcell, job_type_id=16)
    artic_results_path = make_results_directory_artic(
        flowcell.id, artic_task.id, allow_create=False
    )
    if check:
        return artic_results_path.exists()
    if barcodeName:
        # The path to the coverage array for this barcode
        coverage_path = Path(
            f"{artic_results_path}/{barcodeName}/{barcodeName}.coverage.dat"
        )
        return flowcell, artic_results_path, artic_task.id, coverage_path
    return flowcell, artic_results_path, artic_task.id, ""


def remove_duplicate_sequences_numpy(array, master=False, minimum=None):
    """
        Remove stretches of the reference with the same values
        :param array: The array with a value for each reference location
        :param master: Whether this data is for the master chart or not - if it is and
        the array is longer than 50000 elements we want to downsample it until it matches 50000,
        which is the number of points we're happy to draw
        :param minimum: The minimum data value on the x axis
        :return: a tuple of a list with an x,y coordinate for plotting, the max value on the x axis,
        and the max value on the y axis
    """

    # Check we haven't accidentally sent a list
    if not isinstance(array, np.ndarray):
        print(
            f"Expecting Numpy array, received {type(array)}, converting "
            f"to np.ndarray and removing duplicate sequences."
        )
        array = np.array(array)

    # get width of selected window, starting at x coord and finishing at x coord + window width
    if minimum is not None:
        # We aren't starting at x = 0 so add the values on
        xmax = minimum + array.size
        index_start = minimum
    else:
        xmax = array.size
        index_start = 1
    # Get largest value
    ymax = array.max()
    ymin = array.min()
    a = array
    # create an index position, so we can draw the point above the correct base on the graph
    index_position = np.arange(index_start, index_start + a.size, dtype=np.uint32)
    # Remove all the indexes where there is duplicated numbers in sequence
    a = np.insert(a, 0, -1)
    a = np.append(a, -1)
    index_position = index_position[a[2:] != a[:-2]]
    # remove all the duplicate sequential values from the array of actual data values
    a = array[a[2:] != a[:-2]]
    a[a == np.inf] = 0
    a[a == 1.7976931348623157e308] = 0
    # Create a list of tuples, with the x coord as first tuple element and the y value as second
    x_y_values = list(zip(index_position, a))
    # if this is for the coverage master chart we need 50000 values or less for the sake of the browser
    if len(x_y_values) > 50000 and master:
        x_y_values = sampler(x_y_values)
    return x_y_values, xmax, ymax, ymin


def sampler(array):
    """
    If the array is too large > 50000 points, take a sample using steps until we have 50000 or less points
    :param array: A list of tuple values for expected benefit steps
    :type array: list
    :return: a smaller sampled numpy array
    """
    # Get the step size for slicing, getting the array to 50000 points
    step_size = int(len(array) / 50000)
    # Slice is step size
    sampled_array = array[0::step_size]
    return sampled_array


def check_artic_static_exists():
    """
    Check the artic Static directory exists on startup
    Returns
    -------
    None
    """
    if not (Path(BASE_DIR) / "artic" / "static").exists():
        (Path(BASE_DIR) / "artic" / "static").mkdir()
        print("Artic Doesn't exist")
    if not (Path(BASE_DIR) / "artic" / "static" / "artic").exists():
        (Path(BASE_DIR) / "artic" / "static" / "artic").mkdir()


def make_results_directory_artic(flowcell_id, task_id, allow_create=True):
    """
    Make a results directory
    Parameters
    ----------
    flowcell_id: int
        Primary key of the flowcell entry in the database
    task_id: int
        Primary key of the task record in the database.
    allow_create: bool
        Allow the creation of the directory if it doesn't already exist

    Returns
    -------
    results_dir: pathlib.PosixPath
        PosixPath pointing to the results directory
    """
    environmental_results_directory = get_env_variable("MT_ARTIC_RESULTS_DIR")
    artic_dir = Path(f"{environmental_results_directory}/artic/")
    if not artic_dir.exists() and allow_create:
        Path.mkdir(artic_dir)
    results_dir = Path(f"{environmental_results_directory}/artic/Temp_results")
    if not results_dir.exists() and allow_create:
        Path.mkdir(results_dir)
    results_dir = Path(
        f"{environmental_results_directory}/artic/Temp_results/{flowcell_id}_{task_id}_artic"
    )
    if not results_dir.exists() and allow_create:
        Path.mkdir(results_dir)
    return results_dir


def get_amplicon_stats(
    amplicon_band_coords,
    num_amplicons,
    flowcell_id,
    barcode_name,
    time_stamp,
    read_count,
    coverage_array=None,
):
    """
    Get statistics about each amplicon across the coverage array, such as # failed, partial and succesful, mean, std dev
    Parameters
    ----------
    amplicon_band_coords: list of tuple
        list of tuples of start stop across the genome
    num_amplicons: int
        The number of amplicons, equivalent to len(amplicon_band_coords)
    flowcell_id: int
        The primary key of this flowcell
    barcode_name: str
        The name of the barcode
    time_stamp: datetime.datetime
        The time stamp for all these metrics
    read_count: int
        The read count of total mapped reads for this barcode
    coverage_array: np.ndarray
        The array of coverage across the genome

    Returns
    -------
    collection.NamedTuple
    """
    Stats = namedtuple(
        "Stats",
        [
            "time_stamp",
            "coverage",
            "partial_amplicon_count",
            "success_amplicon_count",
            "failed_amplicon_count",
            "read_count",
            "mean_of_amplicon_means",
            "std_dev",
            "variance",
            "would_fire",
            "amplicon_coverage_means",
            "amplicon_coverage_medians",
        ],
    )
    a = np.array(amplicon_band_coords)[:, :2]
    a = a.astype(np.uint32)
    (
        flowcell,
        artic_results_path,
        artic_task_id,
        coverage_path,
    ) = quick_get_artic_results_directory(flowcell_id, barcode_name)
    if coverage_array is None:
        try:
            with open(coverage_path, "rb") as fh:
                coverage = np.fromfile(fh, dtype=np.uint16)
        except FileNotFoundError as e:
            raise e
    else:
        coverage = coverage_array
    amplicon_coverages_median = []
    amplicon_coverages_mean = []
    failed_amplicon_count = 0
    partial_amplicon_count = 0
    cant_count = 0
    for bin_start, bin_end in a:
        amplicon_coverage = coverage[bin_start:bin_end]
        amplicon_coverages_mean.append(np.mean(amplicon_coverage))
        amplicon_median_coverage = np.median(amplicon_coverage)
        amplicon_coverages_median.append(amplicon_median_coverage)
        # For some amplicon schemes we end up with a situation where amplicons basically don't exits so need to chack for nan
        if math.isnan(amplicon_median_coverage):
            cant_count += 1
        elif int(amplicon_median_coverage) == 0:
            failed_amplicon_count += 1
        elif int(amplicon_median_coverage) < 20:
            partial_amplicon_count += 1
    successful_amplicon_counts = (
        num_amplicons - partial_amplicon_count - failed_amplicon_count - cant_count
    )
    ## We need to pick up when we have nan and replace with 0s
    amplicon_mean_array = np.nan_to_num(np.array(amplicon_coverages_mean))
    amplicon_median_array = np.nan_to_num(np.array(amplicon_coverages_median))
    mean_of_amplicon_means = round(amplicon_mean_array.mean(), 2)
    std_dev = round(amplicon_mean_array.std(), 2)
    variance = round(amplicon_mean_array.var(), 2)
    mean_coverage = round(coverage.mean(), 2)
    would_fire = (
        amplicon_mean_array[amplicon_mean_array > 20].size / num_amplicons * 100
    ) > 90
    barcode_stats = Stats(
        time_stamp=time_stamp,
        coverage=mean_coverage,
        partial_amplicon_count=partial_amplicon_count,
        success_amplicon_count=successful_amplicon_counts,
        failed_amplicon_count=failed_amplicon_count,
        read_count=read_count,
        mean_of_amplicon_means=mean_of_amplicon_means,
        std_dev=std_dev,
        variance=variance,
        would_fire=would_fire,
        amplicon_coverage_means=amplicon_mean_array,
        amplicon_coverage_medians=amplicon_median_array,
    )
    return barcode_stats


def predict_barcode_will_finish(
    amplicon_median_array, num_barcodes, total_mapped_reads_count
):
    """
    Predict if after
    Parameters
    ----------
    amplicon_median_array: np.ndarray
        Median coverage for each amplicon
    num_barcodes: int
        The total number of barcodes detected in the run
    total_mapped_reads_count:
        Total number of reads that have mapped to n-cov in this task

    Returns
    -------
    bool
        If this equation thinks the barcode will finish given enough time sequencing

    """
    # 100000 reads per barcode in a run
    ideal_reads_count_constant = int(get_env_variable("MT_IDEAL_READ_CONSTANT"))
    minimum_required_amplicons = int(get_env_variable("MT_ARTIC_MIN_AMPS_PERC"))
    coverage_per_amplicon = int(get_env_variable("MT_COVERAGE_PER_AMPLICON"))
    predicted_coverages = (
        amplicon_median_array
        / total_mapped_reads_count
        * num_barcodes
        * ideal_reads_count_constant
    )
    return (
        predicted_coverages[predicted_coverages > coverage_per_amplicon].size
        / amplicon_median_array.size
    ) * 100 > minimum_required_amplicons


@app.task()
def update_pangolin():
    """

    Returns
    -------

    """
    print("Updating pangolin")
    logger.info("Updating pangolin")
    cmd = [
        "bash",
        "-c",
        f"source {get_env_variable('MT_CONDA_PREFIX')} && conda activate pangolin && pangolin --update",
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    print(out)
    logger.info(out)
    print(repr(err))
    logger.error(err)
