import datetime
import json,gzip
from collections import defaultdict
from urllib.parse import parse_qs

import pandas as pd
from django.shortcuts import render
from rest_framework import status
from rest_framework.decorators import api_view, renderer_classes
from rest_framework.renderers import JSONRenderer, TemplateHTMLRenderer
from rest_framework.response import Response

from alignment.models import PafSummaryCov
# Create your views here.
from artic.models import ArticBarcodeMetadata, ArticFireConditions
from artic.task_artic_alignment import np
from artic.utils import (
    get_amplicon_band_data,
    quick_get_artic_results_directory,
    remove_duplicate_sequences_numpy,
    get_all_results,
    get_amplicon_stats,
)
from minknow_data.models import Flowcell
from minotourapp.utils import get_env_variable
from reads.models import JobMaster, FlowcellSummaryBarcode, Barcode
from reference.models import ReferenceInfo


@api_view(["GET", "POST"])
def fire_conditions_list(request, pk):
    """
    Get the fire conditions for the artic pipeline for this flowcell
    Parameters
    ----------
    request: rest_framework.request.Request
    pk: int
        The primary key of the flowcell record

    Returns
    -------

    """
    if request.method == "GET":
        try:
            data = ArticFireConditions.objects.filter(flowcell_id=pk).values()
            return Response({"data": data}, status=status.HTTP_200_OK)
        except ArticFireConditions.DoesNotExist:
            return Response({"data": []}, status=status.HTTP_404_NOT_FOUND)
    elif request.method == "POST":
        percent_amplicons = request.data.get("num-amplicons-input", False)
        x_coverage = request.data.get("x-coverage-input", False)
        if not all([percent_amplicons, x_coverage]):
            return (
                "Error in request - missing percent amplicons or x coverage choice.",
                status.HTTP_400_BAD_REQUEST,
            )
        ArticFireConditions.objects.create(
            flowcell_id=pk,
            percent_of_amplicons=percent_amplicons,
            x_coverage=x_coverage,
        )
        return Response(
            "Artic firing conditions updated.", status=status.HTTP_201_CREATED
        )


@api_view(["DELETE"])
def fire_condition_detail(request, pk, f_id):
    """
    detail for deleting a specific fire condition by its id
    Parameters
    ----------
    request: rest_framework.request.Request
        The Ajax request body
    pk: int
        The primary key of the flowcell for the artic task
    f_id: int
        Primary key of fire condition to be deleted

    Returns
    -------
    str
        Message string to be displayed

    """
    if request.method == "DELETE":
        ArticFireConditions.objects.get(pk=f_id).delete()
        return Response(
            "Condition was deleted successfully.", status=status.HTTP_204_NO_CONTENT
        )


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
    flowcell, artic_results_path, artic_task_id, _ = quick_get_artic_results_directory(
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
    # Whether to log the data
    # 1 is true 0 is false
    log_coverage = int(request.GET.get("logCoverage", 0))
    if not flowcell_id or not barcode:
        return Response(
            "Flowcell ID and Barcode Name are required.",
            status=status.HTTP_400_BAD_REQUEST,
        )
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
    if log_coverage:
        coverage = coverage.astype(np.float16)
        coverage[coverage == 0] = 0.1
    # Remove duplicate elements in series
    x_y_cov, xmax_cov, ymax_cov, ymin_cov = remove_duplicate_sequences_numpy(
        coverage, True
    )

    x_y_cov = np.array(x_y_cov).round(1)
    # The return dictionary
    data_dict = {
        "coverage": {"xmax": xmax_cov, "ymax": ymax_cov, "data": x_y_cov},
    }
    return Response(data_dict)


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
    log_coverage = int(request.GET.get("logCoverage", 0))
    flowcell_id = request.GET.get("flowcellId")
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
    # We're gonna use the numpy memmap to read only the part of the array we need
    # Create a key matching the contents of this value
    coverage = np.memmap(coverage_path, np.uint16)[mini:maxi]
    if log_coverage:
        coverage = coverage.astype(np.float16)
        coverage[coverage == 0] = 0.1
    # # Remove the duplicate values until a change in the array
    x_y_cov, xmax_cov, ymax_cov, ymin_cov = remove_duplicate_sequences_numpy(
        coverage, minimum=mini
    )
    # Get the primer bands
    scheme = "nCoV-2019"
    scheme_version = "V3"
    amplicon_band_coords, colours = get_amplicon_band_data(scheme, scheme_version)
    data_dict = {
        "coverage": {"xmax": xmax_cov, "ymax": ymax_cov, "data": x_y_cov},
        "amplicon_bands": amplicon_band_coords,
        "colour_lookup": colours,
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
        return Response(
            "no coverage tracking task running on this flowcell.",
            status=status.HTTP_400_BAD_REQUEST,
        )
    try:
        return_data = defaultdict(list)
        queryset = (
            PafSummaryCov.objects.filter(job_master__id=artic_task.id)
            .exclude(barcode_name="unclassified")
            .values_list("barcode_name", "average_read_length", "read_count")
        )
        for barcode, read_length, read_count in queryset:
            return_data["barcodes"].append(barcode)
            return_data["average_read_length"].append(read_length)
            return_data["read_counts"].append(read_count)
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
        return Response(
            "Please specify a flowcellId parameter.",
            status=status.HTTP_400_BAD_REQUEST,
        )
    flowcell, artic_results_path, jm_id, _ = quick_get_artic_results_directory(flowcell_id)
    artic_task = JobMaster.objects.get(pk=jm_id)
    if not artic_task:
        return Response(
            "no coverage tracking task running on this flowcell.",
            status=status.HTTP_400_BAD_REQUEST,
        )
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
            flowcell_id=flowcell_id, job_master=artic_task
        )
    }
    # dictionaries change in place
    scheme = "nCoV-2019"
    scheme_version = "V3"
    amplicon_band_coords, colours = get_amplicon_band_data(scheme, scheme_version)
    num_amplicons = len(amplicon_band_coords)
    time_stamp = datetime.datetime.now()
    for paf_summary_cov in queryset:
        barcode_name = paf_summary_cov["barcode_name"]
        amplicon_stats = get_amplicon_stats(
            amplicon_band_coords,
            num_amplicons,
            flowcell_id,
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
                    artic_results_path / barcode_name / f"{barcode_name}_ARTIC_medaka.csv.gz"
                )
                VoCs=VoCs_df['phe-label']
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
        paf_summary_cov["VoC-Warn"] = VoCs
        amp_stats_dict = amplicon_stats._asdict()
        amp_stats_dict.pop("amplicon_coverage_medians")
        amp_stats_dict.pop("amplicon_coverage_means")
        paf_summary_cov.update(amp_stats_dict)
        paf_summary_cov["lineage"] = lineage
        paf_summary_cov["projected_to_finish"] = artic_metadata[barcode_name].projected_to_finish
    if not queryset:
        return Response(
            f"No coverage summaries found for this task {artic_task.id}.",
            status=status.HTTP_204_NO_CONTENT,
        )
    return Response({"data": queryset})

@api_view(("GET",))
@renderer_classes((TemplateHTMLRenderer, JSONRenderer))
def get_artic_voc_html(request):
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
        return Response(
            "No flowcell ID or barcode provided.", status=status.HTTP_400_BAD_REQUEST
        )
    ## Now we check to see if a json report exists for this flowcell and barcode.
    flowcell, artic_results_path, jm_id, _ = quick_get_artic_results_directory(flowcell_id)
    json_path = artic_results_path / selected_barcode / "json_files" /  f"{selected_barcode}_ARTIC_medaka.json.gz"
    if json_path.exists():
        with gzip.open(json_path, 'r') as fin:
            data = json.loads(fin.read().decode('utf-8'))
            #datastring = json.dumps(data)
            #datastring = datastring.replace("-","_")
            #data = json.loads(datastring)
            data = hyphen_to_underscore(data)
    else:
        data = {}
    csv_path = artic_results_path / selected_barcode / "csv_files" / f"{selected_barcode}_ARTIC_medaka.csv"
    if csv_path.exists():
        df = pd.read_csv(
            csv_path
        )
        html_string = df.T.to_html(classes="table table-sm table-responsive", border=0, justify="left")
        data["hidden_html_string"] = html_string

    csv_path = artic_results_path / selected_barcode / f"{selected_barcode}_ARTIC_medaka.csv.gz"
    if csv_path.exists():
        df = pd.read_csv(
            csv_path
        )
        html_string = df.to_html(classes="table table-sm table-responsive", border=0, index=False, justify="left")
        data["hidden_html_string2"] = html_string

    #vcf_path = artic_results_path / selected_barcode / f"{selected_barcode}.annotated.vcf"

    #if vcf_path.exists():
    #    data["hidden_html_string3"] = "found it"

    lineage_path = artic_results_path / selected_barcode / "lineage_report.csv.gz"
    if lineage_path.exists():
        try:
            data["lineage"] = pd.read_csv(
                artic_results_path / selected_barcode / "lineage_report.csv.gz"
            )["lineage"][0]
        except FileNotFoundError:
            data["lineage"] = "Unknown"


    return render(
        request,
        "artic-variant-of-concern.html",
        context={"artic_barcode_VoC": data},
    )

def hyphen_to_underscore(dictionary):
    """
    Takes an Array or dictionary and replace all the hyphen('-') in any of its keys with a underscore('_')
    :param dictionary:
    :return: the same object with all hyphens replaced by underscore
    """
    # By default return the same object
    final_dict = dictionary

    # for Array perform this method on every object
    if type(dictionary) is type([]):
        final_dict = []
        for item in dictionary:
            final_dict.append(hyphen_to_underscore(item))

    # for dictionary traverse all the keys and replace hyphen with underscore
    elif type(dictionary) is type({}):
        final_dict = {}
        for k, v in dictionary.items():
            # If there is a sub dictionary or an array perform this method of it recursively
            if type(dictionary[k]) is type({}) or type(dictionary[k]) is type([]):
                value = hyphen_to_underscore(v)
                final_dict[k.replace('-', '_')] = value
            else:
                final_dict[k.replace('-', '_')] = v

    return final_dict


@api_view(("GET",))
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
        return Response(
            "No flowcell ID or barcode provided.", status=status.HTTP_400_BAD_REQUEST
        )
    # see if we have a command waiting to be run
    try:
        artic_command_jm = bool(
            JobMaster.objects.get(
                job_type_id=17, barcode__name=selected_barcode, flowcell_id=flowcell_id
            )
        )
    except JobMaster.DoesNotExist:
        artic_command_jm = False
    orm_object = ArticBarcodeMetadata.objects.filter(
        flowcell_id=flowcell_id, barcode__name=selected_barcode
    ).last()
    if not orm_object:
        return Response("No data found", status=status.HTTP_404_NOT_FOUND)
    # First iteration we may may not have FlowcellSumamryBarcodes, so cal
    if not orm_object.percentage_of_reads_in_barcode:
        try:
            barcode_numbers = FlowcellSummaryBarcode.objects.get(
                flowcell_id=flowcell_id, barcode_name=selected_barcode
            )
            all_numbers = FlowcellSummaryBarcode.objects.filter(
                flowcell_id=flowcell_id, barcode_name="All reads"
            ).values_list("read_count", flat=True)
            total_reads = 0
            for all_number in all_numbers:
                total_reads += all_number
            orm_object.percentage_of_reads_in_barcode = (
                barcode_numbers.read_count / total_reads * 100
            )
        except FlowcellSummaryBarcode.DoesNotExist as e:
            orm_object.percentage_of_reads_in_barcode = 0
    # [[new key, old key]]
    new_key_names = [
        ["Avg. Coverage", "average_coverage"],
        ["Var. Coverage", "variance_coverage"],
        ["Min. Coverage", "minimum_coverage"],
        ["Max. Coverage", "maximum_coverage"],
        ["% reads in run", "percentage_of_reads_in_barcode"],
        ["Has Finished", "has_finished"],
        ["Has Sufficient Coverage", "has_sufficient_coverage"],
    ]
    results_files = [
        ("Consensus sequence", "consensus"),
        ("Box plot", "box-plot"),
        ("Bar plot", "bar-plot"),
        ("Fail VCF", "fail-vcf"),
        ("Pass VCF", "pass-vcf"),
        ("Input fasta", "input-fasta"),
        ("Pangolin lineages CSV", "pangolin-lineages"),
        ("Sorted Bam", "sorted-bam"),
        ("Sorted Bam index", "sorted-bam-bai"),
    ]
    old_dict = orm_object.__dict__
    context_dict = {key[0]: old_dict[key[1]] for key in new_key_names}
    context_dict["hidden_barcode_pk"] = orm_object.barcode.id
    context_dict["hidden_barcode_name"] = orm_object.barcode.name
    context_dict["hidden_flowcell_id"] = flowcell_id
    context_dict["hidden_job_master_id"] = orm_object.job_master.id
    context_dict["hidden_results_files"] = results_files
    context_dict["hidden_has_finished"] = old_dict["has_finished"]
    context_dict["hidden_has_suff"] = old_dict["has_sufficient_coverage"]
    context_dict["hidden_marked_for_rerun"] = old_dict["marked_for_rerun"]
    context_dict["hidden_destroy_evidence"] = bool(
        int(get_env_variable("MT_DESTROY_ARTIC_EVIDENCE"))
    )
    context_dict["hidden_triggered_by_cleanup"] = (
        orm_object.has_finished and not orm_object.has_sufficient_coverage
    )
    context_dict["hidden_has_command_job_master"] = artic_command_jm
    (
        flowcell,
        artic_results_path,
        artic_task_id,
        _,
    ) = quick_get_artic_results_directory(flowcell_id)
    fastq_path = artic_results_path / selected_barcode / f"{selected_barcode}.fastq"
    fastq_path_gz = fastq_path.with_suffix(".fastq.gz")
    context_dict["hidden_has_fastq"] = fastq_path.exists() or fastq_path_gz.exists()
    if context_dict["hidden_has_finished"]:
        csv_path =artic_results_path / selected_barcode / "lineage_report.csv.gz"
        if csv_path.exists():
            df = pd.read_csv(
                artic_results_path / selected_barcode / "lineage_report.csv.gz"
            )
            html_string = df.T.to_html(classes="table-responsive table-striped", border=0)
            context_dict["hidden_html_string"] = html_string
    return render(
        request,
        "artic-barcode-metadata.html",
        context={"artic_barcode_metadata": context_dict},
    )


@api_view(["POST"])
def manually_create_artic_command_job_master(request):
    """
    Manually create a Job Master for the run_artic_command for a given barcode.
    Parameters
    ----------
    request: rest_framework.request.Request
        Django rest framework request object. Should contain flowcellId, job_master_id, barcode_name as data.
    Returns
    -------
    status: int
        The HTTP code of the request.
    """
    job_type_id = request.data.get("jobTypeId", None)
    barcode_pk = request.data.get("barcodePk", None)
    flowcell_id = request.data.get("flowcellId", None)
    job_master_pk = request.data.get("jobPk", None)
    if not job_type_id or not barcode_pk or not flowcell_id:
        return Response(
            "Flowcell id, barcode id, jobTypeId are required fields.",
            status=status.HTTP_400_BAD_REQUEST,
        )
    try:
        reference = JobMaster.objects.get(pk=job_master_pk).reference
        barcode = Barcode.objects.get(pk=barcode_pk)

    except [ReferenceInfo.DoesNotExist, ReferenceInfo.MultipleObjectsReturned] as e:
        return Response(
            f"Exception: {e}", status=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )
    # TODO potential bug here where the barcode name on the reads is not the barcode name we save

    job_master, created = JobMaster.objects.get_or_create(
        job_type_id=job_type_id,
        reference=reference,
        barcode=barcode,
        flowcell_id=flowcell_id,
    )
    if created:
        return Response(
            "Successfully created task to be consumed by Celery later.",
            status=status.HTTP_200_OK,
        )
    else:
        return Response(
            "Already a JobMaster for this barcode on this flowcell.",
            status=status.HTTP_403_FORBIDDEN,
        )


@api_view(["GET"])
def get_results_package(request):
    """
    TODO clear results after download
    Create a results package for artic barcode data, for user download. Need a way to clear it. Maybe a timed call? A celery task?
    Parameters
    ----------
    request: rest_framework.request.Request
        Post request data. Contains files chosen, flowcell id and barcode.
    Returns
    -------
        Status

    """
    params = json.loads(request.GET.get("params", None))
    all_barcodes = params.get("all", False)
    flowcell_id = params.get("flowcellId", None)
    barcode_pk = int(params.get("selectedBarcode", None))
    chosen = parse_qs(request.GET.get("string"))
    barcode_name = Barcode.objects.get(pk=barcode_pk).name if barcode_pk != -1 else None
    barcode_name = "" if all_barcodes else barcode_name
    if not flowcell_id or not barcode_name and not all_barcodes:
        return Response(
            "Flowcell Id or barcode required in request.",
            status=status.HTTP_400_BAD_REQUEST,
        )

    try:
        flowcell = Flowcell.objects.get(pk=flowcell_id)
    except Flowcell.DoesNotExist:
        return Response("Flowcell not found.", status=status.HTTP_404_NOT_FOUND)

    # TODO add check to see if task is complete
    _, path_to_results, _, _ = quick_get_artic_results_directory(
        flowcell_id, barcode_name
    )
    response = get_all_results(path_to_results, flowcell, barcode_name, chosen.keys())
    return response


@api_view(["GET"])
@renderer_classes((TemplateHTMLRenderer, JSONRenderer))
def get_results_modal_html(request, pk):
    """
    return the html for the modal for the all results download functionality
    Parameters
    ----------
    request: rest_framework.request.Request
        The ajax request body
    pk: int
        The primary key of the flowcell object in the database

    Returns
    -------
    html

    """
    results_files_extra = [
        ("Input fasta", "input-fasta"),
        ("Sorted Bam", "sorted-bam"),
        ("Sorted Bam Index", "sorted-bam-bai"),
    ]
    results_files = [
        ("Consensus sequence", "consensus"),
        ("Box plot", "box-plot"),
        ("Bar plot", "bar-plot"),
        ("Fail VCF", "fail-vcf"),
        ("Pass VCF", "pass-vcf"),
        ("Pangolin lineages CSV", "pangolin-lineages"),
    ]
    if not int(get_env_variable("MT_DESTROY_ARTIC_EVIDENCE")):
        results_files.extend(results_files_extra)
    context_dict = {"hidden_results_files": results_files}
    return render(request, "all-results-modal.html", context={"context": context_dict})


@api_view(["GET"])
@renderer_classes((TemplateHTMLRenderer, JSONRenderer))
def png_html(request):
    """
    Return the html for the PNG cards, and give them the context to render the images
    Parameters
    ----------
    request: rest_framework.request.Request
        Request object, contains flowcell primary key and barcode for sample to be shown

    Returns
    -------

    """
    flowcell_id = request.GET.get("flowcellId", None)
    selected_barcode = request.GET.get("selectedBarcode", None)
    if not flowcell_id or not selected_barcode:
        return Response(
            "No flowcell ID or barcode provided.", status=status.HTTP_400_BAD_REQUEST
        )
    orm_object = ArticBarcodeMetadata.objects.filter(
        flowcell_id=flowcell_id, barcode__name=selected_barcode
    ).last()
    if not orm_object:
        return Response("No data found", status=status.HTTP_404_NOT_FOUND)
    if not orm_object.has_finished:
        return Response("No results yet.", status=status.HTTP_204_NO_CONTENT)
    (
        flowcell,
        artic_results_path,
        artic_task_id,
        _,
    ) = quick_get_artic_results_directory(flowcell_id)
    context_list = [
        [
            # TODO we need to add the flowcell_pk_job_master_pk folder before here
            f"artic/{flowcell_id}_{artic_task_id}_artic/{selected_barcode}-barplot.png",
            f"{selected_barcode} Artic Bar Plot",
            "This plot is generated directly by the artic pipeline which has been run on this data set. It illustrates the coverage of each amplicon for this barcode (using pass reads only) at the time the pipeline was run.",
        ],
        [
            f"artic/{flowcell_id}_{artic_task_id}_artic/{selected_barcode}-boxplot.png",
            f"{selected_barcode} Artic Box Plot",
            "Distribution of amplicon coverage in the amplicon pools at the time this pipeline was run. This plot is automatically generated by the artic pipline.",
        ],
    ]
    png_dict = {"srcs": context_list, "has_finished": orm_object.has_finished}
    return render(request, "artic-pngs.html", context={"png": png_dict})


@api_view(["GET"])
def get_amplicon_bands_for_master(request):
    """
    TODO should contain scheme and scheme version
    Get the amplicon x and y coordinates for the master chart on initialisation
    Parameters
    ----------
    request: rest_framework.request.Request
        The request object.
    Returns
    -------
    list of lists
        List of the start and stop coordinates, and colour lookup

    """
    # Get the primer bands
    scheme = "nCoV-2019"
    scheme_version = "V3"
    amplicon_band_coords, colours = get_amplicon_band_data(scheme, scheme_version)
    data_dict = {"amplicon_band_coords": amplicon_band_coords, "colours": colours}
    return Response(data_dict, status=status.HTTP_200_OK)


@api_view(["PATCH"])
def rerun_artic_command(request):
    """
    Manually set a completed fire artic command jobmaster to be rerun
    Parameters
    ----------
    request: rest_framework.request.Request
        The request object
    Returns
    -------

    """
    flowcell_id = request.data.get("flowcellId", None)
    selected_barcode_pk = request.data.get("selectedBarcodePk", None)
    if not flowcell_id or not selected_barcode_pk:
        return Response(
            "No flowcell ID or barcode provided.", status=status.HTTP_400_BAD_REQUEST
        )
    orm_object = ArticBarcodeMetadata.objects.filter(
        flowcell_id=flowcell_id, barcode_id=selected_barcode_pk
    ).last()
    artic_command_task = JobMaster.objects.filter(
        flowcell_id=flowcell_id, job_type_id=17, barcode_id=selected_barcode_pk
    ).last()
    if not artic_command_task.complete:
        return Response("Command not yet run.", status=status.HTTP_403_FORBIDDEN)
    artic_command_task.complete = False
    artic_command_task.save()
    orm_object.marked_for_rerun = True
    orm_object.save()
    return Response(
        "Task sucessfully listed for re-run", status=status.HTTP_204_NO_CONTENT
    )


@api_view(["GET"])
def get_artic_pie_chart_data(request):
    """
    Get pie chart data for the artic page
    Parameters
    ----------
    request: rest_framework.request.Request
        Rest framework request object. Contains Flowcell ID
    Returns
    -------

    """
    flowcell_id = request.GET.get("flowcellId", None)
    if not flowcell_id:
        return Response("No flowcell ID provided.", status=status.HTTP_400_BAD_REQUEST)
    artic_task = JobMaster.objects.filter(
        flowcell_id=flowcell_id, job_type__name="Track Artic Coverage"
    ).last()
    queryset = PafSummaryCov.objects.filter(job_master__id=artic_task.id).values_list(
        "barcode_name", "read_count"
    )
    if not queryset:
        return Response("No data.", status=status.HTTP_204_NO_CONTENT)
    total = 0
    classified_total = 0
    unclassified_total = 0
    for barcode, read_count in queryset:
        total += read_count
        if barcode == "unclassified":
            unclassified_total += read_count
        else:
            classified_total += read_count

    classified_total = round(classified_total / total * 100, 2)
    unclassified_total = round(unclassified_total / total * 100, 2)
    series_data = [
        {
            "name": "Proportion of Reads Classified",
            "ColorByPoint": True,
            "data": [
                {
                    "name": "Classified",
                    "y": classified_total,
                    "sliced": True,
                    "selected": True,
                },
                {"name": "Unclassified", "y": unclassified_total},
            ],
        }
    ]
    return Response(series_data, status=status.HTTP_200_OK)


@api_view(["PATCH"])
def run_all_incomplete(request, pk):
    """
    Manually run all incomplete barcodes on an artic task through the artic pipeline
    Parameters
    ----------
    request: rest_framework.request.Request
        The Ajax request object
    pk: int
        Flowcell primary key

    Returns
    -------
    str
        Message pertaining to success or failure of view
    """
    jm = JobMaster.objects.get(flowcell_id=pk, job_type_id=16)
    incomplete_barcodes = jm.JobMastersArticBarcodeMetadatas.filter(has_finished=False)
    for barcode in incomplete_barcodes:
        job_master, created = JobMaster.objects.get_or_create(
            job_type_id=17,
            reference=jm.reference,
            barcode=barcode.barcode,
            flowcell=jm.flowcell,
        )
    return Response("Marked all incomplete barcodes for rerun.")


@api_view(["PATCH"])
def mark_all_barcodes_for_pipeline(request, pk):
    """
    Mark all barcodes under an artic task to be run through the pipeline,
     irrespective of whether they have been run before or not.
    Parameters
    ----------
    request: rest_framework.request.Request
        The Ajax request object
    pk: int
        The primary key of the object.

    Returns
    -------
    str
        Message pertaining to success or failure of the view.
    """
    jm = JobMaster.objects.get(flowcell_id=pk, job_type_id=16)
    barcodes = jm.JobMastersArticBarcodeMetadatas.all()
    for barcode in barcodes:
        job_master, created = JobMaster.objects.get_or_create(
            job_type_id=17,
            reference=jm.reference,
            barcode=barcode.barcode,
            flowcell=jm.flowcell,
        )
        if not created:
            job_master.complete = False
            job_master.save()
            barcode.marked_for_rerun = True
            barcode.save()
    return Response("Marked all barcodes to be run.")
