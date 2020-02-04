import json
from collections import defaultdict

import numpy as np
import pandas as pd
from django.db.models import Sum
from django.http import HttpResponse, JsonResponse
from rest_framework.decorators import api_view
from alignment.models import PafRoughCov, PafSummaryCov
from alignment.utils import calculate_coverage_new
from rest_framework.response import Response
from rest_framework import status

from reads.models import JobMaster


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):

        if isinstance(obj, np.integer):
            return int(obj)

        elif isinstance(obj, np.floating):
            return float(obj)

        else:
            return super(MyEncoder, self).default(obj)


@api_view(["GET"])
def rough_coverage_complete_chromosome_flowcell(
    request, task_id, barcode_name, read_type_id, chromosome_id
):
    """
    Fetch data for coverage master charts
    Parameters
    ----------
    request: rest_framework.request.Request
        Django rest framework request object.
    task_id: int
        The PK of the job master database row.
    barcode_name: str
        The name of the barcode.
    read_type_id: int
        The PK of the read type database entry.
    chromosome_id: int
        The PK of the chromosome.

    Returns
    -------

    """
    print(task_id, barcode_name, read_type_id, chromosome_id, request.user)

    queryset = (
        PafRoughCov.objects.filter(job_master__id=task_id)
        .filter(flowcell__owner=request.user)
        .filter(barcode_name=barcode_name)
        .filter(chromosome__id=chromosome_id)
        .filter(read_type__id=read_type_id)
        .values("p")
        .annotate(Sum("i"))
        .order_by("p")
    )

    print(queryset)
    print(queryset)

    queryset_df = pd.DataFrame.from_records(queryset.values("p", "i__sum"))
    print(queryset_df)

    queryset_df["i__sum__cumsum"] = queryset_df["i__sum"].cumsum()

    result = queryset_df[["p", "i__sum__cumsum"]]

    # this method return the json between double quotes and requires a JSON.parse on the client side
    return JsonResponse(result.to_json(orient="values"), safe=False)


@api_view(["GET"])
def flowcell_paf_alignment_list(
    request, task_id, barcode_name, read_type_id, chromosome_id, start, end
):
    # def flowcell_paf_alignment_list(request, flowcell_id, barcodegroup_id, read_type_id, chromosome_id, start, end):

    min_extreme = request.GET.get("min", "")
    max_extreme = request.GET.get("max", "")

    result_list = calculate_coverage_new(
        task_id, barcode_name, read_type_id, chromosome_id, start, end
    )

    return JsonResponse(json.dumps(result_list, cls=NumpyEncoder), safe=False)


@api_view(["GET"])
def paf_summary_json(request, pk):
    """
    Return the data for the minimap2 coverage summary table
    Parameters
    ----------
    request: rest_framework.Request.request
        The request object that the table sends with the parameters for this data.
    pk: int
        The primary key of the flowcell.

    Returns:
    -------
    Response: rest_framework.Response.response
        The Response object, containing response data, status code.
    result: dict
        The resultant data, contain the values for each column in the DataTables table.

    """

    query_columns = [
        "barcode_name",
        "chromosome__line_name",
        "reference_line_length",
        "read_count",
        "total_length",
        "average_read_length",
        "coverage",
    ]

    draw = int(request.GET.get("draw", 0))

    start = int(request.GET.get("start", 0))

    length = int(request.GET.get("length", 10))

    end = start + length
    # Which column s
    search_value = request.GET.get("search[value]", "")

    order_column = request.GET.get("order[0][column]", "")
    # ascending descending
    order_dir = request.GET.get("order[0][dir]", "")

    if not search_value == "":

        queryset = PafSummaryCov.objects.filter(job_master__flowcell_id=pk).filter(
            chromosome__line_name__icontains=search_value
        )

    else:

        queryset = PafSummaryCov.objects.filter(job_master__flowcell_id=pk)

    print(queryset)
    print(queryset.values())

    records_total = queryset.count()

    queryset = queryset[start:end]

    df = pd.DataFrame.from_records(
        queryset.values(
            "barcode_name",
            "chromosome__line_name",

            "job_master__id",
            "read_count",
            "total_length",
            "reference_line_length",
        )
    )

    df["average_read_length"] = (
        df["total_length"].div(df["read_count"]).round(decimals=3)
    )

    df["coverage"] = (
        df["total_length"].div(df["reference_line_length"]).round(decimals=3)
    )

    if order_column:

        if order_dir == "desc":
            df = df.sort_values(query_columns[int(order_column)], ascending=False)

        else:
            df = df.sort_values(query_columns[int(order_column)], ascending=True)

    dictdf = df.to_dict("records")

    result = {
        "draw": draw,
        "recordsTotal": records_total,
        "recordsFiltered": records_total,
        "data": dictdf,
    }

    return Response(result, status=200)


@api_view(["GET"])
def flowcell_paf_summary_cov(request, pk):
    """
    Return the data used to draw the charts on the mapping summary page.
    Parameters
    ----------
    request: rest_framework.request.Request
        The rest framework Request object.
    pk:int
        The primary key of the flowcell this data is for.

    Returns
    -------
    dict
        Contains the data for the charts.
    """
    queryset = PafSummaryCov.objects.filter(job_master__flowcell_id=pk)

    response = []

    for record in queryset:

        paf_summary_cov = {
            "id": record.job_master.id,
            "barcode_name": record.barcode_name,
            "read_type_name": "Template",
            "reference_line_name": record.chromosome.line_name,
            "read_count": record.read_count,
            "total_length": record.total_length,
            "chrom_cover": round(record.total_length / record.reference_line_length, 2),
            "avg_read_len": round(record.total_length / record.read_count),
        }

        response.append(paf_summary_cov)

    return JsonResponse({"data": response})


@api_view(["GET"])
def flowcellreferences_used_by_run(request, flowcell_id):
    """
    Get the references that have been used in any mapping jobs on a flowcell,
    return them to populate the dropdown box.
    Parameters
    ----------
    request: rest_framework.request.Request
        Django rest framework request object.
    flowcell_id: int
        The PK of the flowcell database row.

    Returns
    -------
    json
        The data for the drop down.
    """
    references = (
        PafSummaryCov.objects.filter(job_master__flowcell__id=flowcell_id)
        .values(
            "job_master__id",
            "job_master__reference__name",
            "job_master__reference__id",
            "barcode_name",
            "chromosome__line_name",
            "chromosome__id",
            "read_type__id",
            "read_type__name"
        )
        .distinct()
    )

    result = [
        {
            "task_id": r["job_master__id"],
            "flowcell_id": flowcell_id,
            "chromosome_name": r["chromosome__line_name"],
            "chromosome_id": r["chromosome__id"],
            "barcode_name": r["barcode_name"],
            "reference_name": r["job_master__reference__name"],
            "reference_id": r["job_master__reference__id"],
            "read_type_name": r["read_type__name"],
            "read_type_id": r["read_type__id"]
        }
        for r in references
    ]

    return HttpResponse(json.dumps(list(result)), content_type="application/json")


@api_view(["GET"])
def get_coverage_summary(request, pk):
    """
    Get a list of all the chromosomes and their coverage.
    Parameters
    ----------
    request: rest_framework.request.Request
        Django rest framework object
    pk: int
        The primary key of the flowcell
    Returns
    -------
    If names is True - Dict of Reference to Contig, Contig to Barcodes, else return coverage per contig

    """
    print(request.GET)
    # if the request is from the Dropdown for notifications
    if request.GET.get("names", False):
        queryset = JobMaster.objects.filter(flowcell__id=pk, job_type__name="minimap2").values_list(
            "reference__name", "reference__referencelines__line_name", "flowcell__runs__barcodes__name",
            "reference__id", "reference__referencelines__id", "flowcell__runs__barcodes__id")

        if not queryset:
            return Response("No mapping tasks available on this flowcell.", status=400)

        reference_to_contig_dict = defaultdict(list)

        # create a dict of string keys and values to create the conditional drop downs
        # contig[0] - reference name
        # contig[1] - contig names
        # contig[2] - barcode names
        # contig[3] - reference pk
        # contig[4] - contig id
        # contig[5] - barcodes id
        for contig in queryset:

            contig_to_barcode_dict = defaultdict(list)

            reference_to_contig_dict[contig[0]] = [contig_to_barcode_dict, contig[3]]

            if not contig_to_barcode_dict[contig[1]]:

                contig_to_barcode_dict[contig[1]].append(contig[4])

            contig_to_barcode_dict[contig[1]].append((contig[2], contig[5]))

        return Response(reference_to_contig_dict, status=200)

    queryset = PafSummaryCov.objects.filter(job_master__flowcell__id=pk)

    if not queryset:
        return Response(
            "No alignment results found for this flowcell",
            status=status.HTTP_204_NO_CONTENT,
        )

    df = pd.DataFrame.from_records(
        queryset.values(
            "barcode_name",
            "chromosome__line_name",
            "total_length",
            "reference_line_length",
        )
    )

    df["coverage"] = (
        df["total_length"].div(df["reference_line_length"]).round(decimals=3)
    )

    df = df.drop(columns=["reference_line_length", "total_length"])

    df.set_index("chromosome__line_name", inplace=True)

    dictdf = df.to_dict("index")

    result = {"data": dictdf}

    return Response(result, status=200)
