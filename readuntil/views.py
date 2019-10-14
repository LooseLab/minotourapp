"""
Read Until views that do the work
"""
from django.db.models import Q
from rest_framework.decorators import api_view
from rest_framework.response import Response
from pathlib import Path
from jobs.models import JobMaster
from reads.models import Flowcell, Run
import numpy as np
from readuntil.utils import get_or_create_results_directory
from readuntil.models import (
    ExpectedBenefitChromosomes,
    RejectedFastqRead,
    AcceptedFastqRead,
)
import pandas as pd


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
    index_position = np.arange(
        index_start, index_start + a.size, dtype=np.uint32
    )
    # Remove all the indexes where there is duplicated numbers in sequence

    a = np.insert(a, 0, -1)

    a = np.append(a, -1)

    index_position = index_position[a[2:] != a[:-2]]
    # remove all the duplicate sequential values from the array of actual data values
    a = array[a[2:] != a[:-2]]

    a[a == np.inf] = 0

    a[a == 1.7976931348623157e+308] = 0

    # Create a list of tuples, with the x coord as first tuple element and the y value as second
    x_y_values = list(zip(index_position, a))

    # if this is for the coverage master chart we need 50000 values or less for the sake of the browser
    if len(x_y_values) > 50000 and master:
        print("Down sampling for master")
        x_y_values = sampler(x_y_values)

    return x_y_values, xmax, ymax, ymin


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


@api_view(["GET"])
def get_benefit_data_master(request):
    """
    Return data for mock master plot, and the yAxis maximum for the detail charts
    :param request: rest framework http request object
    :type request: rest_framework.request.Request
    :return:
    """
    flowcell_id = request.GET.get("flowcellId")

    latest_eb_task = (
        JobMaster.objects.filter(
            job_type__name="ExpectedBenefit", flowcell_id=flowcell_id
        )
        .last()
        .id
    )

    basepath = get_or_create_results_directory(flowcell_id, latest_eb_task)

    chromosome = request.GET.get("chromosome_chosen", "NC_003210.1")

    # The path to the coverage array for this chromosome
    coverage_path = Path(
        f"{basepath}/coverage_{chromosome}_{flowcell_id}_{latest_eb_task}.dat"
    )

    try:
        with open(coverage_path, "rb") as fh:
            coverage = np.fromfile(fh, dtype=np.uint16)

    except FileNotFoundError as e:
        Response("File not found", status=404)
        raise e
    # Remove duplicate elements in series
    x_y_cov, xmax_cov, ymax_cov, ymin_cov = remove_duplicate_sequences_numpy(
        coverage, True
    )

    # The return dictionary
    data_dict = {
        "coverage": {"xmax": xmax_cov, "ymax": ymax_cov, "data": x_y_cov}
    }

    return Response(data_dict)


@api_view(["GET"])
def get_benefit_data_detail(request):
    """
    The complete benefit data for the detail chart from the selected regions of the master chart
    :param request: Contains the min and max x axis values in the body
    :type request: rest_framework.request.Request
    :return: a dictionary of the data from the pickle for all the detail charts
    """
    # TODO Refactor to be in for loop so it's DRYER
    # The min and max values of the x axis
    # the minimum and maxiumum x coordinates for the detail graph
    mini = int(request.GET.get("min"))
    maxi = int(request.GET.get("max"))

    flowcell_id = request.GET.get("flowcellId")

    print("REQUEST {}".format(request.GET))
    #print (request.GET.get("chromosome_chosen"))

    chromosome = request.GET.get("chromosome_chosen", "NC_003210.1")

    jm = JobMaster.objects.filter(
        job_type__name="ExpectedBenefit", flowcell_id=flowcell_id
    ).last()

    if jm is None:
        return Response("Expected Benefit task not found", status=404)

    latest_eb_task = (
        jm.id
    )
    # the base path to the directory that we have stored the binary files in
    basepath = get_or_create_results_directory(flowcell_id, latest_eb_task)
    # the paths to each array
    coverage_path = Path(
        f"{basepath}/coverage_{chromosome}_{flowcell_id}_{latest_eb_task}.dat"
    )

    counts_path = Path(
        f"{basepath}/counts_{chromosome}_{flowcell_id}_{latest_eb_task}.dat"
    )

    benefits_path = Path(
        f"{basepath}/benefits_{chromosome}_{flowcell_id}_{latest_eb_task}.dat"
    )

    mask_path_forward = Path(
        f"{basepath}/mask_forward_{chromosome}_{flowcell_id}_{latest_eb_task}.dat"
    )

    mask_path_reverse = Path(
        f"{basepath}/mask_reverse_{chromosome}_{flowcell_id}_{latest_eb_task}.dat"
    )

    scores_forward_path = Path(
        f"{basepath}/scores_forward_{chromosome}_{flowcell_id}_{latest_eb_task}.dat"
    )

    scores_reverse_path = Path(
        f"{basepath}/scores_reverse_{chromosome}_{flowcell_id}_{latest_eb_task}.dat"
    )

    cost_forward_path = Path(
        f"{basepath}/cost_forward_{chromosome}_{flowcell_id}_{latest_eb_task}.dat"
    )

    cost_reverse_path = Path(
        f"{basepath}/cost_reverse_{chromosome}_{flowcell_id}_{latest_eb_task}.dat"
    )

    fixed_ben_forward_path = Path(
        f"{basepath}/fixed_benefits_forward{chromosome}_{flowcell_id}_{latest_eb_task}.dat"
    )

    fixed_ben_reverse_path = Path(
        f"{basepath}/fixed_benefits_reverse{chromosome}_{flowcell_id}_{latest_eb_task}.dat"
    )

    # the file path and data type we need as a tuple
    file_path_list = [
        (coverage_path, np.uint16),
        (
            counts_path,
            [
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
        ),
        (benefits_path, np.float64),
        (mask_path_forward, np.int8),
        (mask_path_reverse, np.int8),
        (scores_forward_path, np.float64),
        (scores_reverse_path, np.float64),
        (cost_forward_path, np.float64),
        (cost_reverse_path, np.float64),
        (fixed_ben_forward_path, np.float64),
        (fixed_ben_reverse_path, np.float64),
    ]
    # use these strings as a key
    file_contents = [
        "coverage",
        "counts",
        "benefits",
        "mask_forward",
        "mask_reverse",
        "scores_forward",
        "scores_reverse",
        "costs_forward",
        "costs_reverse",
        "fixed_ben_forward",
        "fixed_ben_reverse",
    ]
    # Get a list with a dictionary of all the numpy arrays
    data_dict = {}

    # We're gonna use the numpy memmap to read only the part of the array we need
    for index, file_path_or_dtype in enumerate(file_path_list):
        # Create a key matching the contents of this value
        data_dict[file_contents[index]] = np.memmap(
            file_path_or_dtype[0], dtype=file_path_or_dtype[1]
        )[mini:maxi]
        print()

    # # Remove the duplicate values until a change in the array
    x_y_cov, xmax_cov, ymax_cov, ymin_cov = remove_duplicate_sequences_numpy(
        data_dict["coverage"], minimum=mini
    )
    # # Remove the duplicate values until a change in the array, on array with numpy nans removed
    x_y_benefit, xmax_benefit, ymax_benefit, ymin_benefit = remove_duplicate_sequences_numpy(
        np.nan_to_num(data_dict["benefits"], copy=False), minimum=mini
    )
    # # Remove the duplicate values until a change in the array
    x_y_fwd_mask, xmax_fwd_mask, ymax_fwd_mask, ymin_fwd_mask = remove_duplicate_sequences_numpy(
        data_dict["mask_forward"], minimum=mini
    )
    print(xmax_cov)
    # Remove the duplicate values until a change in the array
    x_y_rev_mask, xmax_rev_mask, ymax_rev_mask, ymin_fwd_mask = remove_duplicate_sequences_numpy(
        data_dict["mask_reverse"], minimum=mini
    )

    x_y_fwd_scores, xmax_fwd_scores, ymax_fwd_scores, ymin_fwd_scores = remove_duplicate_sequences_numpy(
        np.nan_to_num(data_dict["scores_forward"], posinf=1), minimum=mini
    )

    x_y_rev_scores, xmax_rev_scores, ymax_rev_scores, ymin_rev_scores = remove_duplicate_sequences_numpy(
        np.nan_to_num(data_dict["scores_reverse"], posinf=1), minimum=mini
    )

    x_y_fwd_costs, xmax_fwd_costs, ymax_fwd_costs, ymin_fwd_costs = remove_duplicate_sequences_numpy(
        np.nan_to_num(data_dict["costs_forward"], copy=False), minimum=mini
    )

    x_y_rev_costs, xmax_rev_costs, ymax_rev_costs, ymin_rev_costs = remove_duplicate_sequences_numpy(
        np.nan_to_num(data_dict["costs_reverse"], copy=False), minimum=mini
    )

    x_y_fwd_fix_ben, xmax_fwd_fix_ben, ymax_fwd_fix_ben, ymin_fwd_fix_ben = remove_duplicate_sequences_numpy(
        np.nan_to_num(data_dict["fixed_ben_forward"], copy=False), minimum=mini
    )

    x_y_rev_fix_ben, xmax_rev_fix_ben, ymax_rev_fix_ben, ymin_fwd_fix_ben = remove_duplicate_sequences_numpy(
        np.nan_to_num(data_dict["fixed_ben_reverse"], copy=False), minimum=mini
    )

    data_dict = {
        "coverage": {"xmax": xmax_cov, "ymax": ymax_cov, "data": x_y_cov},
        "forwardMask": {
            "xmax": xmax_fwd_mask,
            "ymax": ymax_fwd_mask,
            "data": x_y_fwd_mask,
        },
        "reverseMask": {
            "xmax": xmax_rev_mask,
            "ymax": ymax_rev_mask,
            "data": x_y_rev_mask,
        },
        "benefits": {
            "xmax": xmax_benefit,
            "ymax": ymax_benefit,
            "data": x_y_benefit,
        },
        "scoresFwd": {
            "xmax": xmax_fwd_scores,
            "ymax": ymax_fwd_scores,
            "ymin": ymin_fwd_scores,
            "data": x_y_fwd_scores,
        },
        "scoresRev": {
            "xmax": xmax_rev_scores,
            "ymax": ymax_rev_scores,
            "ymin": ymin_rev_scores,
            "data": x_y_rev_scores,
        },
        "costsFwd": {
            "xmax": xmax_fwd_costs,
            "ymax": ymax_fwd_costs,
            "ymin": ymin_fwd_costs,
            "data": x_y_fwd_costs,
        },
        "costsRev": {
            "xmax": xmax_rev_costs,
            "ymax": ymax_rev_costs,
            "ymin": ymin_rev_costs,
            "data": x_y_rev_costs,
        },
        "fixBenFwd": {
            "xmax": xmax_fwd_fix_ben,
            "ymax": ymax_fwd_fix_ben,
            "data": x_y_fwd_fix_ben,
        },
        "fixBenRev": {
            "xmax": xmax_rev_fix_ben,
            "ymax": ymax_rev_fix_ben,
            "data": x_y_rev_fix_ben,
        },
    }

    return Response(data_dict)


@api_view(["GET"])
def get_chromosomes(request):
    """
    Get a list of the chromosomes that we have EB results for
    :param request: The DRF request object from the web client
    :type request: rest_framework.request.Request
    :return:
    """
    flowcell_id = request.GET.get("flowcellId")
    # get the django ORM object for the latest EB task
    latest_eb_task = JobMaster.objects.filter(
        job_type__name="ExpectedBenefit", flowcell_id=flowcell_id
    ).last()
    # Get all the chromosome names we have attached to this task that we have results for
    chromosome_names = ExpectedBenefitChromosomes.objects.filter(
        task=latest_eb_task
    ).values_list("chromosome__line_name", flat=True)

    return Response(list(chromosome_names))


@api_view(["GET", "POST"])
def get_or_create_rejected_reads(request):
    """
    Get or create the read ids for rejected fastq reads in a read until experiment
    :param request: Django rest framework HTTP request object
    :type request: rest_framework.request.Request
    :return: Response of all reject fastq reads for a certain flowcell
    """

    # Lookup for whether we are fetching Accepted or rejected data

    if request.method is "GET":
        # get the flowcell id we are using in this request, sent by the browser
        flowcell_id = request.query_params.get("flowcell_id", -1)

        # no flowcell Id was provided
        if flowcell_id == -1:
            return Response(
                "No flowcell id was provided. Please provide a flowcell ID.",
                status=404,
            )

        flowcell = Flowcell.objects.get(
            Q(name=flowcell_id) | Q(pk=flowcell_id)
        )
        # no flowcell found
        if not flowcell_id:
            return Response(
                "No flowcell was found mathcing that name or ID", status=404
            )

        run_id = request.query_params.get("run_id", -1)

        if run_id == -1:
            return Response(
                "No run_id was provided. Please provide a flowcell ID.",
                status=404,
            )

        run = Run.objects.get(runid=run_id)

        if not run_id:
            return Response("No run was found matching that ID", status=404)

        rejected_read_ids = RejectedFastqRead.objects.filter(
            run=run
        ).values_list("read_id", flat=True)

        return Response(rejected_read_ids, status=200)
    # it's a post so we're getting data from the readUntil client
    else:

        flowcell_id = request.DATA["flowcell_id"]

        run_id = request.DATA["run_id"]

        flowcell = Flowcell.objects.get(
            Q(name=flowcell_id) | Q(pk=flowcell_id)
        )

        if not flowcell:
            return Response(
                "No flowcell was found matching that name or ID", status=404
            )

        run = Run.objects.get(runid=run_id)

        if not run:
            return Response("No run was found matching that ID", status=404)

        # if it's a list, we are saving both accepted and rejected read ids
        fastq_reads = request.DATA["read_ids"]

        if not fastq_reads:
            return Response("No fastq read ids provided", status=500)
        # Create a dataframe of the new read_ids
        df = pd.DataFrame(fastq_reads, columns=["read_id"])
        # broadcast the flowcell objects down the dataframe
        df["flowcell"] = flowcell
        # broadcast the run objects down the dataframe
        df["run"] = run
        # transform to list of dictionaries row wise, keyed by column names
        saving_dict = df.to_dict(orient="records")

        del df

        # bulk create them
        reads_saved = RejectedFastqRead.objects.bulk_create(
            [RejectedFastqRead(**item) for item in saving_dict],
            batch_size=2000,
        )


@api_view(["GET", "POST"])
def get_or_create_accepted_reads(request):
    """
    Get or create the read ids for Accepted fastq reads in a read until experiment
    :param request: Django rest framework HTTP request object
    :type request: rest_framework.request.Request
    :return:
    """
    # Lookup for whether we are fetching Accepted or rejected data

    if request.method is "GET":
        # get the flowcell id we are using in this request, sent by the browser
        flowcell_id = request.query_params.get("flowcell_id", -1)

        # no flowcell Id was provided
        if flowcell_id == -1:
            return Response(
                "No flowcell id was provided. Please provide a flowcell ID.",
                status=404,
            )

        flowcell = Flowcell.objects.get(
            Q(name=flowcell_id) | Q(pk=flowcell_id)
        )
        # no flowcell found
        if not flowcell_id:
            return Response(
                "No flowcell was found mathcing that name or ID", status=404
            )

        run_id = request.query_params.get("run_id", -1)

        if run_id == -1:
            return Response(
                "No run_id was provided. Please provide a flowcell ID.",
                status=404,
            )

        run = Run.objects.get(runid=run_id)

        if not run_id:
            return Response("No run was found matching that ID", status=404)

        accepted_read_ids = AcceptedFastqRead.objects.filter(
            run=run
        ).values_list("read_id", flat=True)

        return Response(accepted_read_ids, status=200)
    # it's a post so we're getting data from the readUntil client
    else:

        flowcell_id = request.DATA["flowcell_id"]

        run_id = request.DATA["run_id"]

        flowcell = Flowcell.objects.get(
            Q(name=flowcell_id) | Q(pk=flowcell_id)
        )

        if not flowcell:
            return Response(
                "No flowcell was found matching that name or ID", status=404
            )

        run = Run.objects.get(runid=run_id)

        if not run:
            return Response("No run was found matching that ID", status=404)

        # if it's a list, we are saving both accepted and rejected read ids
        fastq_reads = request.DATA["read_ids"]

        if not fastq_reads:
            return Response("No fastq read ids provided", status=500)
        # Create a dataframe of the new read_ids
        df = pd.DataFrame(fastq_reads, columns=["read_id"])
        # broadcast the flowcell objects down the dataframe
        df["flowcell"] = flowcell
        # broadcast the run objects down the dataframe
        df["run"] = run
        # transform to list of dictionaries row wise, keyed by column names
        saving_dict = df.to_dict(orient="records")

        del df

        # bulk create them
        reads_saved = AcceptedFastqRead.objects.bulk_create(
            [AcceptedFastqRead(**item) for item in saving_dict],
            batch_size=2000,
        )