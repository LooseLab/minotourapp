"""
Useful functions that don't belong in views or tasks
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from django.shortcuts import get_object_or_404

from artic.task_artic_alignment import make_results_directory_artic
from minotourapp.utils import get_env_variable
from reads.models import Flowcell, JobMaster

colour_palette = ["#ffb3ba", "#ffdfba", "#ffffba", "#baffc9", "#bae1ff"]


def get_amplicon_band_data(scheme, scheme_version):
    """
    Retrieve coordinates on reference for amplicon bands, and a colour scheme for any amplicon pools
    Parameters
    ----------
    scheme: str
        Artic primer scheme to get data from
    scheme_version: str
        The version, if any
    Returns
    -------
    (list of string, dict)
        A list of coordinates, start and stop on x axis, and a colour scheme lookup dictionary for amplicon pools.

    """
    json_file_path = get_amplicon_json_file_path(scheme, scheme_version)
    with open(json_file_path, "r") as fh:
        amplicon_bands = json.load(fh)
    # Get data
    amplicon_band_coords = json.loads(amplicon_bands["amplicons"])
    colours = {
        amplicon_bands["pools"][index]: colour_palette[index]
        for index in range(len(amplicon_bands["pools"]))
    }
    return amplicon_band_coords, colours


def get_amplicon_json_file_path(scheme, scheme_version):
    """
    Get the file path to the amplicon primer JSON file
    Parameters
    ----------
    scheme: str
        Artic primer scheme to get data from
    scheme_version: str
        The version, if any
    Returns
    -------
    pathlib.PosixPath
        File path to the JSON primer file.

    """
    scheme_dir = Path(get_env_variable("MT_ARTIC_SCEHEME_DIR"))
    artic_results_dir_primers = (
        Path(get_env_variable("MT_ARTIC_RESULTS_DIR")) / "artic/" / "primers_files"
    )
    if not (Path(get_env_variable("MT_ARTIC_RESULTS_DIR")) / "artic/").exists():
        (Path(get_env_variable("MT_ARTIC_RESULTS_DIR")) / "artic/").mkdir()
    # Needs to be dynamic here
    bed_file = "nCoV-2019.scheme.bed"
    json_file = f"{scheme}.primers.json"
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
    filepath: pathlib.PosixPath
        Absolute file path to the primer scheme bed file
    json_file: str
        Name to give the JSON file created from the bed file
    artic_results_primer_dir: pathlib.PosixPath
        Absolute path to the directory to store amplicon JSON in

    Returns
    -------
    pathlib.PosixPath
        File path to JSON file, stored in Artic Results Dir
    """
    # TODO not dynamic, how make dynamic, very hardcoded
    # TODO Dropdown on task start for scheme dir?
    print(filepath)
    df = pd.read_csv(filepath, sep="\t", header=None)
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
    json_data = {
        "amplicons": df[["primer_start", "primer_end", 4]].to_json(orient="values"),
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

    artic_task = get_object_or_404(
        JobMaster, flowcell=flowcell, job_type__name="Track Artic Coverage"
    )

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

    return flowcell, artic_results_path, artic_task.id


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
        print("Down sampling for master")
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