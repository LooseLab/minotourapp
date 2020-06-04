from collections import defaultdict
from pathlib import Path

from django.shortcuts import get_object_or_404, render
from rest_framework import status
from rest_framework.decorators import api_view, renderer_classes
from rest_framework.renderers import JSONRenderer, TemplateHTMLRenderer
from rest_framework.response import Response

from alignment.models import PafSummaryCov
# Create your views here.
from artic.models import ArticBarcodeMetadata
from artic.task_artic_alignment import make_results_directory_artic, np
from reads.models import Flowcell, JobMaster, FlowcellSummaryBarcode


def quick_get_artic_results_directory(flowcell_id, barcodeName="", check=False):
    """

    Parameters
    ----------
    flowcell_id : int
        Primary key of flowcell record that we are looking
    barcodeName : str optional
        The barcode that the data is from. If provided coverage path is also returned.
    check : bool
        Check whether we have results or not.

    Returns
    -------
    (flowcell: reads.models.Flowcell, artic_results_path: pathlib.PosixPath, artic_task.id : int)



    """
    flowcell = Flowcell.objects.get(pk=flowcell_id)

    artic_task = get_object_or_404(JobMaster,
        flowcell=flowcell, job_type__name="Track Artic Coverage"
    )

    artic_results_path = make_results_directory_artic(
        flowcell.id, artic_task.id, allow_create=False
    )

    if check:
        return artic_results_path.exists()

    if barcodeName:
        # TODO Note this chromsome is hard coded here because we are only doing this on Covid at the moment
        chromosome = "NC_045512.2"
        # The path to the coverage array for this barcode
        coverage_path = Path(
            f"{artic_results_path}/{barcodeName}/coverage_{chromosome}_{flowcell_id}_{artic_task.id}_{barcodeName}.dat"
        )
        return flowcell, artic_results_path, artic_task.id, coverage_path

    return flowcell, artic_results_path, artic_task.id


@api_view(["GET"])
def get_artic_barcodes(request):
    """
    Return a list of Barcodes that have Artic coverage associated with them.
    Parameters
    ----------
    request: rest_framework.request.Request

    Returns
    -------

    """
    flowcell_id = int(request.GET.get("flowcellId", None))

    if not flowcell_id:
        return Response(
            "Please provide a flowcell id.", status=status.HTTP_400_BAD_REQUEST
        )

    flowcell, artic_results_path, artic_task_id = quick_get_artic_results_directory(
        flowcell_id
    )

    barcodes_list = []

    if not artic_results_path.exists():
        return Response(
            "Artic results directory not found.", status=status.HTTP_404_NOT_FOUND
        )

    for item in artic_results_path.iterdir():
        if item.is_dir():
            barcodes_list.append(item.stem)

    if not barcodes_list:
        return Response("No results found.", status=status.HTTP_404_NOT_FOUND)

    return Response(barcodes_list, status=status.HTTP_200_OK)


@api_view(["GET"])
def get_artic_master_chart_data(request):
    """
    Get the master chart data for the selected barcode of a flowcell
    Parameters
    ----------
    request

    Returns
    -------

    """
    flowcell_id = request.GET.get("flowcellId", None)

    barcode = request.GET.get("barcodeChosen", None)

    if not flowcell_id or not barcode:
        return Response("Flowcell ID and Barcode Name are required.", status=status.HTTP_400_BAD_REQUEST)

    (
        flowcell,
        artic_results_path,
        artic_task_id,
        coverage_path,
    ) = quick_get_artic_results_directory(flowcell_id, barcode)

    try:
        with open(coverage_path, "rb") as fh:
            coverage = np.fromfile(fh, dtype=np.uint16)

    except FileNotFoundError as e:
        Response("Coverage file not found.", status=status.HTTP_404_NOT_FOUND)
        raise e
    # Remove duplicate elements in series
    x_y_cov, xmax_cov, ymax_cov, ymin_cov = remove_duplicate_sequences_numpy(
        coverage, True
    )

    # The return dictionary
    data_dict = {"coverage": {"xmax": xmax_cov, "ymax": ymax_cov, "data": x_y_cov}}

    return Response(data_dict)


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


@api_view(["GET"])
def get_artic_detail_chart_data(request):
    """
    The complete benefit data for the detail chart from the selected regions of the master chart

    Parameters
    ----------
    request: rest_framework.request.Request
        The request object recieved by Django Rest Framework

    Returns
    -------

    """
    # The min and max values of the x axis
    # the minimum and maxiumum x coordinates for the detail graph
    mini = int(request.GET.get("min"))
    maxi = int(request.GET.get("max"))

    flowcell_id = request.GET.get("flowcellId")

    print("REQUEST {}".format(request.GET))
    # print (request.GET.get("chromosome_chosen"))

    barcode = request.GET.get("barcodeChosen", "NOTFOUND")

    (
        flowcell,
        artic_results_path,
        artic_task_id,
        coverage_path,
    ) = quick_get_artic_results_directory(flowcell_id, barcode)

    if artic_task_id is None:
        return Response(
            "Artic task not found for this flowcell.", status=status.HTTP_404_NOT_FOUND
        )

    if not coverage_path.exists():
        return Response(
            "No results files found for this flowcell.",
            status=status.HTTP_404_NOT_FOUND,
        )

    # the file path and data type we need as a tuple
    file_path_list = [
        (coverage_path, np.uint16),
    ]
    # use these strings as a key
    file_contents = [
        "coverage",
    ]
    # Get a list with a dictionary of all the numpy arrays
    data_dict = {}

    # We're gonna use the numpy memmap to read only the part of the array we need
    for index, file_path_or_dtype in enumerate(file_path_list):
        # Create a key matching the contents of this value
        data_dict[file_contents[index]] = np.memmap(
            file_path_or_dtype[0], dtype=file_path_or_dtype[1]
        )[mini:maxi]

    # # Remove the duplicate values until a change in the array
    x_y_cov, xmax_cov, ymax_cov, ymin_cov = remove_duplicate_sequences_numpy(
        data_dict["coverage"], minimum=mini
    )

    data_dict = {
        "coverage": {"xmax": xmax_cov, "ymax": ymax_cov, "data": x_y_cov},
    }

    return Response(data_dict)


@api_view(["GET"])
def get_artic_column_chart_data(request):
    """
    Return the data for the Average Coverage per barcode and Average Read length per barcode
    column charts for a specified flowcell.
    Parameters
    ----------
    request: rest_framework.request.Request
        Django rest framework request class. Contains Flowcell ID as a param..
    Returns
    -------

    """

    flowcell_id = request.GET.get("flowcellId", None)

    if not flowcell_id:
        return Response(
            "No flowcellId specified. Please check request.",
            status=status.HTTP_400_BAD_REQUEST,
        )
    artic_task = JobMaster.objects.filter(
        flowcell_id=flowcell_id, job_type__name="Track Artic Coverage"
    ).last()

    if not artic_task:
            return Response("no coverage tracking task running on this flowcell.", status=status.HTTP_400_BAD_REQUEST)

    try:
        return_data = defaultdict(list)
        queryset = PafSummaryCov.objects.filter(job_master__id=227).values_list("barcode_name", "average_read_length", "coverage")
        for barcode, read_length, coverage in queryset:
            return_data["barcodes"].append(barcode)
            return_data["average_read_length"].append(read_length)
            return_data["coverage"].append(coverage)

    except TypeError as e:
        print(e.args)
        return Response(e, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    if not queryset:
        return Response(
            "No values for this flowcell found. Please check task has been created.",
            status=status.HTTP_404_NOT_FOUND,
        )

    return Response(return_data, status=status.HTTP_200_OK)


@api_view(["GET"])
def get_artic_summary_table_data(request):
    """
    Return the PafSummaryCov data for a given flowcell and populate a datatable
    Parameters
    ----------
    request: rest_framework.request.Request
        Request object that contains the flowcell Primary key
    Returns
    -------

    """
    flowcell_id = request.GET.get("flowcellId", None)

    if not flowcell_id:
        return Response("Please specify a flowcell_id parameter.", status=status.HTTP_400_BAD_REQUEST)

    artic_task = JobMaster.objects.filter(
        flowcell_id=flowcell_id, job_type__name="Track Artic Coverage"
    ).last()

    if not artic_task:
        return Response("no coverage tracking task running on this flowcell.", status=status.HTTP_400_BAD_REQUEST)

    queryset = PafSummaryCov.objects.filter(job_master=artic_task).values("barcode_name",
                                                                          "chromosome__line_name",
                                                                          "average_read_length",
                                                                          "coverage",
                                                                          "job_master__id",
                                                                          "read_count",
                                                                          "total_length",
                                                                          "reference_line_length", )

    if not queryset:
        return Response(f"No coverage summaries found for this task {artic_task.id}.", status=status.HTTP_204_NO_CONTENT)

    return Response({"data": queryset})


@api_view(('GET',))
@renderer_classes((TemplateHTMLRenderer, JSONRenderer))
def get_artic_barcode_metadata_html(request):
    """

    Parameters
    ----------
    request: rest_framework.request.Request
        Request body, params: the flowcell PK and selected barcode

    Returns
    -------

    """
    flowcell_id = request.GET.get("flowcellId", None)
    selected_barcode = request.GET.get("selectedBarcode", None)
    if not flowcell_id or not selected_barcode:
        return Response("No flowcell ID or barcode provided.", status=status.HTTP_400_BAD_REQUEST)
    orm_object = ArticBarcodeMetadata.objects.filter(flowcell_id=flowcell_id, barcode__name=selected_barcode).last()

    if not orm_object:
        return Response("No data found", status=status.HTTP_404_NOT_FOUND)

    # First iteration we may may not have FlowcellSumamryBarcodes
    if not orm_object.percentage_of_reads_in_barcode:
        try:
            barcode_numbers = FlowcellSummaryBarcode.objects.get(flowcell_id=flowcell_id, barcode_name=selected_barcode)
            all_numbers = FlowcellSummaryBarcode.objects.filter(flowcell_id=flowcell_id, barcode_name="All reads").values_list("read_count", flat=True)
            total_reads = 0
            for all_number in all_numbers:
                total_reads += all_number
            orm_object.percentage_of_reads_in_barcode = barcode_numbers.read_count / total_reads

        except FlowcellSummaryBarcode.DoesNotExist as e:
            orm_object.percentage_of_reads_in_barcode = 0

            # [[new key, old key]]
    new_key_names = [["Avg. Coverage", "average_coverage"], ["Var. Coverage", "variance_coverage"], ["Min. Coverage", "minimum_coverage"], ["Max. Coverage", "maximum_coverage"],
                     ["% in barcode", "percentage_of_reads_in_barcode"]]
    old_dict = orm_object.__dict__
    context_dict = {key[0]: old_dict[key[1]] for key in new_key_names}
    context_dict["barcode_name"] = orm_object.barcode.name

    try:
        context_dict["has_run"] = JobMaster.objects.get(flowcell_id=flowcell_id, job_type__name="Run Artic",
                                                        barcode__name=selected_barcode).complete
        context_dict["has_reached_sufficient_coverage"] = True
    except JobMaster.DoesNotExist:
        context_dict["has_run"] = False
        context_dict["has_reached_sufficient_coverage"] = False


    return render(request, "artic-barcode-metadata.html", context={"artic_barcode_metadata": context_dict})
