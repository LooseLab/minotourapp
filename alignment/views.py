import json

import numpy as np
import pandas as pd
from django.db.models import F, Sum, Avg
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from alignment.models import PafRoughCov, PafSummaryCov


@api_view(["GET"])
def rough_coverage_complete_chromosome_flowcell(
    request, task_id, barcode_id, read_type_id, chromosome_id
):
    """
    Fetch data for coverage master charts
    Parameters
    ----------
    request: rest_framework.request.Request
        Django rest framework request object.
    task_id: int
        The PK of the job master database row.
    barcode_id: int
        The primary key of the barcode.
    read_type_id: int
        The PK of the read type database entry.
    chromosome_id: int
        The PK of the chromosome.

    Returns
    -------

    """
    queryset = (
        PafRoughCov.objects.filter(job_master__id=task_id)
        .filter(flowcell__owner=request.user)
        .filter(barcode_id=barcode_id)
        .filter(chromosome__id=chromosome_id)
        .filter(read_type__id=read_type_id)
        .order_by("bin_position_end")
    ).values_list("bin_position_start", "bin_coverage")
    # TODO limits to just one reference here when fetching length, could be an issue in the future
    length = (
        PafSummaryCov.objects.filter(job_master_id=task_id, chromosome_id=chromosome_id)
        .first()
        .chromosome.chromosome_length
    )
    return Response(
        {"chartData": queryset, "refLength": length}, status=status.HTTP_200_OK
    )


@api_view(["GET"])
def paf_summary_table_json(request, pk):
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
    queryset = PafSummaryCov.objects.filter(job_master__flowcell_id=pk).exclude(
        job_master__job_type_id=16
    )
    df = pd.DataFrame.from_records(
        queryset.values(
            "barcode_name",
            "chromosome__line_name",
            "chromosome__reference__name",
            "job_master_id",
            "read_count",
            "total_yield",
            "reference_line_length",
        )
    )
    df["average_read_length"] = (
        df["total_yield"].div(df["read_count"]).round(decimals=3)
    )
    df["coverage"] = (
        df["total_yield"].div(df["reference_line_length"]).round(decimals=3)
    )
    dictdf = df.to_dict("records")
    result = {
        "data": dictdf,
    }
    return Response(result, status=200)


@api_view(["GET"])
def flowcell_paf_summary_cov(request, pk):
    """
    Return the data used to draw the per barcode column charts on the mapping summary page.
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
    chromosome_pk = request.GET.get("chromosomeId", None)
    queryset = PafSummaryCov.objects.filter(
        job_master__flowcell_id=pk, chromosome_id=chromosome_pk
    ).exclude(job_master__job_type_id=16)
    if not queryset:
        return Response(
            "No data currently held in database.", status=status.HTTP_404_NOT_FOUND
        )
    df = pd.DataFrame.from_records(
        queryset.values(
            "job_master_id",
            "barcode_name",
            "chromosome_id",
            "chromosome__line_name",
            "reference_line_length",
            "read_count",
            "total_yield",
            "read_type_id",
            "read_type__name",
            "coverage",
            "average_read_length",
        )
    )
    # Weird group by iterations I wrote in a fugue to make the categories right for grouped columns,
    # inserting correct number of 0s at start for each reference
    df["coverage"] = df["coverage"].apply(lambda x: [x])
    df["average_read_length"] = df["average_read_length"].apply(lambda x: [x])
    for index, name in enumerate(df["barcode_name"].unique()):
        if index > 0:
            df.loc[df["barcode_name"] == name, "coverage"] = df[
                df["barcode_name"] == name
            ]["coverage"].apply(lambda x: ([0] * index) + x)
            df.loc[df["barcode_name"] == name, "average_read_length"] = df[
                df["barcode_name"] == name
            ]["average_read_length"].apply(lambda x: ([0] * index) + x)

    df["coverage_series"] = df[
        [
            "chromosome__line_name",
            "coverage",
            "average_read_length",
            "barcode_name",
            "read_type__name",
        ]
    ].T.apply(
        lambda col: {"name": col[3] + " " + col[4] + " " + col[0], "data": col[1]}
    )
    df["average_read_length_series"] = df[
        [
            "chromosome__line_name",
            "coverage",
            "average_read_length",
            "barcode_name",
            "read_type__name",
        ]
    ].T.apply(
        lambda col: {"name": col[3] + " " + col[4] + " " + col[0], "data": col[2]}
    )
    categories = np.unique(df["chromosome__line_name"].values)
    return Response(
        {
            "coverageData": df["coverage_series"].values,
            "avgRLData": df["average_read_length_series"].values,
            "categories": categories,
        }
    )


@api_view(["GET"])
def mapped_references_by_flowcell_list(request, flowcell_id):
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
    if not flowcell_id:
        return Response(
            "No flowcell ID slug provided", status=status.HTTP_400_BAD_REQUEST
        )
    references = (
        PafSummaryCov.objects.filter(job_master__flowcell__id=flowcell_id)
        .exclude(job_master__job_type_id=16)
        .values(
            jmId=F("job_master_id"),
            referenceName=F("job_master__reference__name"),
            referenceId=F("job_master__reference_id"),
            barcodeName=F("barcode_name"),
            barcodeId=F("barcode_id"),
            chromosomeName=F("chromosome__line_name"),
            chromosomeId=F("chromosome_id"),
            readTypeId=F("read_type_id"),
            readTypeName=F("read_type__name"),
        )
    )
    if not references:
        return Response("No data found", status=status.HTTP_204_NO_CONTENT)
    return Response(references, status=status.HTTP_200_OK)


@api_view(["GET"])
def per_genome_coverage_summary(request, flowcell_pk):
    """
    Get the per genome (i.e chromosome level) coverage fot the top column charts
    Parameters
    ----------
    request: rest_framework.request.Request
        The request rest framework object
    flowcell_pk: int
        The flowcell primary key in the database

    Returns
    -------
    list of dict

    """
    queryset = PafSummaryCov.objects.filter(job_master__flowcell_id=flowcell_pk).values("chromosome__line_name").annotate(
        Sum("total_yield"), Avg("average_read_length"))
    chromosome_coverage = []
    chromosome_average_read_length = []
    for chromosome in queryset:
        chromosome_coverage.append(
            {"name": chromosome["chromosome__line_name"], "data": [chromosome["total_yield__sum"]]})
        chromosome_average_read_length.append(
            {"name": chromosome["chromosome__line_name"], "data": [chromosome["average_read_length__avg"]]})
    results = {"coverageData": chromosome_coverage, "avgRLData": chromosome_average_read_length}
    return Response(results, status=status.HTTP_200_OK)
