"""
Read Until views that do the work
"""
from rest_framework.decorators import api_view
from rest_framework.response import Response
from pathlib import Path
from jobs.models import JobMaster
import numpy as np
from readuntil.utils import make_results_directory
from readuntil.models import ExpectedBenefitChromosomes


def find_files_of_type(file_or_directory, file_extensions):
    """Return a list of pathlib.Path of files with chosen extensions
    Parameters
    ----------
    file_or_directory : str
        filepath or a directory
    file_extensions : list
        A list of lowercase file extensions including '.' e.g. ['.txt', '.paf']
    Returns
    -------
    list
        If files with extension are found return list of pathlib.Path objects, otherwise return empty list
    """
    file_or_directory = Path(file_or_directory).expanduser()
    if (
        file_or_directory.is_file()
        and "".join(file_or_directory.suffixes).lower() in file_extensions
    ):
        return [file_or_directory]
    elif file_or_directory.is_dir():
        return [
            x
            for x in file_or_directory.iterdir()
            if "".join(x.suffixes).lower() in file_extensions
        ]
    else:
        return []


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
    # the last element, which we need to append onto our array with deleted duplicates
    last_element = array[-1]

    a = array
    # we're going to cut this off below, so we have to keep it to add it back on
    last_index = a.size
    # create an index position, so we can draw the point above the correct base on the graph
    index_position = np.arange(index_start, index_start + a.size, dtype=np.uint32)
    # Remove all the indexes where there is duplicated numbers in sequence
    index_position = index_position[:-1][a[1:] != a[:-1]]
    # Add the removed last element back on
    index_position = np.append(index_position, last_index)
    # remove all the duplicate sequential values from the array of actual data values
    a = a[:-1][a[1:] != a[:-1]]
    # Add the truncated last element back on
    a = np.append(a, last_element)
    # Create a list of tuples, with the x coord as first tuple element and the y value as second
    x_y_values = list(zip(index_position, a))

    # if this is for the coverage master chart we need 50000 values or less for the sake of the browser
    if len(x_y_values) > 50000 and master:
        print("Down sampling for master")
        x_y_values = sampler(x_y_values)

    return x_y_values, xmax, ymax


def get_y_axis_max(array):
    """
    Get the maximum value in the array
    :param array: The array of data
    :type array: np.nparray
    :return:
    """
    return array.max()


def sampler(array):
    """
    If the array is too large > 50000 points, take a sample using steps
    :param array: Numpy array
    :return: a smaller sampled numpy array
    """
    # Get the step size for slicing, getting the array to 50000 points
    step_size = int(len(array) / 50000)
    # Slice is step size
    sampled_array = array[0::step_size]
    return sampled_array


@api_view(["GET"])
def get_benefit_data_master(request):
    """
    Return data for mock master plot, and the yAxis maximum for the detail charts
    :param request:
    :return:
    """
    flowcell_id = request.GET.get("flowcellId")

    latest_eb_task = JobMaster.objects.filter(job_type__name="ExpectedBenefit", flowcell_id=flowcell_id).last().id

    basepath = make_results_directory(flowcell_id, latest_eb_task)

    chromosome = request.GET.get("chromosome", "NC_003210.1")

    # The path to the coverage array for this chromosome
    coverage_path = Path(f"{basepath}/coverage_{chromosome}_{flowcell_id}_{latest_eb_task}.dat")

    try:
        with open(coverage_path, "rb") as fh:
            coverage = np.fromfile(fh, dtype=np.uint16)

    except FileNotFoundError as e:
        Response("File not found", status=404)
        raise e
    # Remove duplicate elements in series
    x_y_cov, xmax_cov, ymax_cov = remove_duplicate_sequences_numpy(coverage, True)

    # The return dictionary
    data_dict = {"coverage": {"xmax": xmax_cov, "ymax": ymax_cov, "data": x_y_cov},
                 }

    return Response(data_dict)


@api_view(["GET"])
def get_benefit_data_detail(request):
    """
    The complete benefit data for the detail chart from the selected regions of the master chart
    :param request: Contains the min and max x axis values in the body
    :return: a dictionary of the data from the pickle for all the detail charts
    """
    # TODO just for one chromsome at the moment
    # The min and max values of the x axis
    # the minimum and maxiumum x coordinates for the detail graph
    mini = int(request.GET.get("min"))
    maxi = int(request.GET.get("max"))

    flowcell_id = request.GET.get("flowcellId")

    chromosome = request.GET.get("chromosome", "NC_003210.1")

    latest_eb_task = JobMaster.objects.filter(job_type__name="ExpectedBenefit", flowcell_id=flowcell_id).last().id
    # the base path to the directory that we have stored the binary files in
    basepath = make_results_directory(flowcell_id, latest_eb_task)
    # the paths to each array
    coverage_path = Path(
        f"{basepath}/coverage_{chromosome}_{flowcell_id}_{latest_eb_task}.dat")

    counts_path = Path(
        f"{basepath}/counts_{chromosome}_{flowcell_id}_{latest_eb_task}.dat")

    benefits_path = Path(
        f"{basepath}/benefits_{chromosome}_{flowcell_id}_{latest_eb_task}.dat")

    mask_path_forward = Path(
        f"{basepath}/mask_forward_{chromosome}_{flowcell_id}_{latest_eb_task}.dat")

    mask_path_reverse = Path(
        f"{basepath}/mask_reverse_{chromosome}_{flowcell_id}_{latest_eb_task}.dat")

    # the file path and data type we need as a tuple
    file_path_list = [(coverage_path, np.uint16), (counts_path, [("A", np.uint16), ("C", np.uint16),
                                                                 ("G", np.uint16), ("T", np.uint16), ("D", np.uint16,),
                                                                 ("I", np.uint16), ("IC", np.uint16),
                                                                 ("M", np.uint16), ("U", np.bool)]),
                      (benefits_path, np.float64),
                      (mask_path_forward, np.bool),
                      (mask_path_reverse, np.bool)]
    # use these strings as a key
    file_contents = ["coverage", "counts", "benefits", "mask_forward", "mask_reverse"]
    # Get a list with a dictionary of all the numpy arrays
    data_dict = {}

    # We're gonna use the numpy memmap to read only the part of the array we need
    for index, file_path_or_dtype in enumerate(file_path_list):
        # Create a key matching the contents of this value
        data_dict[file_contents[index]] = np.memmap(file_path_or_dtype[0],
                                                    dtype=file_path_or_dtype[1])[mini:maxi]
        print()

    # # Remove the duplicate values until a change in the array
    x_y_cov, xmax_cov, ymax_cov = remove_duplicate_sequences_numpy(data_dict["coverage"],
                                                                   minimum=mini)
    # # Remove the duplicate values until a change in the array, on array with numpy nans removed
    x_y_benefit, xmax_benefit, ymax_benefit = remove_duplicate_sequences_numpy(
        np.nan_to_num(data_dict["benefits"], copy=False),
        minimum=mini)
    # # Remove the duplicate values until a change in the array
    x_y_fwd_mask, xmax_fwd_mask, ymax_fwd_mask = remove_duplicate_sequences_numpy(data_dict["mask_forward"],
                                                                                  minimum=mini)
    print(xmax_cov)
    # Remove the duplicate values until a change in the array
    x_y_rev_mask, xmax_rev_mask, ymax_rev_mask = remove_duplicate_sequences_numpy(data_dict["mask_reverse"],
                                                                                  minimum=mini)

    data_dict = {"coverage": {"xmax": xmax_cov, "ymax": ymax_cov, "data": x_y_cov},

                 "forwardMask": {"xmax": xmax_fwd_mask, "ymax": ymax_fwd_mask,
                                 "data": x_y_fwd_mask},
                 "reverseMask": {"xmax": xmax_rev_mask, "ymax": ymax_rev_mask,
                                 "data": x_y_rev_mask},
                 "benefits": {"xmax": xmax_benefit, "ymax": ymax_benefit,
                              "data": x_y_benefit}
                 }

    return Response(data_dict)


@api_view(["GET"])
def get_chromosomes(request):
    """
    Get a list of the chromosome we have EB results for
    :param request: The DRF request object from the web client
    :return:
    """
    flowcell_id = request.GET.get("flowcellId")
    # get the django ORM object for the latest EB task
    latest_eb_task = JobMaster.objects.filter(job_type__name="ExpectedBenefit", flowcell_id=flowcell_id).last()
    # Get all the chromosome names we have attached to this task that we have results for
    chromosome_names = ExpectedBenefitChromosomes.objects.filter(
        task=latest_eb_task).values_list(
        "chromosome__line_name", flat=True)

    return Response(list(chromosome_names))



