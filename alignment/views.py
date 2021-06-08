import json

import numpy as np
import pandas as pd
from django.db.models import F, Sum, Avg
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from alignment.models import PafSummaryCov, MattsAmazingAlignmentSum
from reads.models import Barcode


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
    test = MattsAmazingAlignmentSum.objects.filter(job_master__id=task_id,
            flowcell__owner=request.user,
            barcode_id=barcode_id,
            chromosome_pk=chromosome_id,
            read_type_id=read_type_id,)
    testDF = pd.DataFrame.from_records(test.values())
    func = lambda s: np.fromstring(s, dtype=int, sep=",")
    testDF['Numpy_bin_position_start'] = testDF['bin_position_start_str'].str.strip("[]").apply(func)
    testDF['Numpy_bin_change'] = testDF['bin_coverage_str'].str.strip("[]").apply(func)
    # testDF['Numpy_bin_position_start'] = testDF['bin_position_start_str'].apply(lambda x: np.array(eval(x)))
    # testDF['Numpy_bin_change'] = testDF['bin_coverage_str'].apply(lambda x: np.array(eval(x)))
    testarray = testDF[['Numpy_bin_position_start', 'Numpy_bin_change', 'rejected_barcode_id']]
    testarray = testarray.set_index(['rejected_barcode_id']).apply(pd.Series.explode).reset_index().astype('int64')
    ###This will generate the old style data we need for now:
    testarray2 = testarray.groupby(['Numpy_bin_position_start'], as_index=False).agg(
                {'Numpy_bin_change': np.sum})
    testarray2['cumsum'] = testarray2['Numpy_bin_change'].transform(np.cumsum)
    testarray2 = testarray2.rename(columns={'Numpy_bin_position_start': 'bin_position_start', 'cumsum': 'bin_coverage'})


    ### This generates the new data set which will have a labelled dictionary.
    testarray = testarray.groupby(['rejected_barcode_id', 'Numpy_bin_position_start'], as_index=False).agg(
        {'Numpy_bin_change': np.sum})
    testarray['cumsum'] = testarray.groupby(['rejected_barcode_id'])['Numpy_bin_change'].transform(np.cumsum)
    testarray = testarray.drop(columns=['Numpy_bin_change'])
    testarray = testarray.set_index('rejected_barcode_id')
    testarray = testarray.rename(columns={'Numpy_bin_position_start': 'bin_position_start', 'cumsum': 'bin_coverage'})

    results_dict = dict()
    for name, group in testarray.groupby(['rejected_barcode_id']):
        mybarcode = Barcode.objects.get(pk=name)
        results_dict[mybarcode.name] = group.values.tolist()

    queryset = testarray2[['bin_position_start', 'bin_coverage']].values
    new_data = results_dict
    sum_to_check=queryset.sum()

    # TODO limits to just one reference here when fetching length, could be an issue
    #  in the future displaying multiple rferences on a plot
    length = (
        PafSummaryCov.objects.filter(job_master_id=task_id, chromosome_pk=chromosome_id)
        .first()
        .reference_line_length
    )
    return Response(
        {"newChartData": new_data, "chartData": queryset, "refLength": length, "sumToCheck": sum_to_check},
        status=status.HTTP_200_OK,
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
    queryset = (
        PafSummaryCov.objects.filter(job_master__flowcell_id=pk)
        .exclude(job_master__job_type_id=16)
        .values(
            "barcode_name",
            "chromosome_name",
            "reference_name",
            "job_master_id",
            "read_count",
            "total_yield",
            "reference_line_length",
            "coverage",
            "average_read_length",
        )
    )
    result = {
        "data": queryset,
    }
    return Response(result, status=200)


@api_view(["GET"])
def per_barcode_coverage_summary(request, pk):
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
    barcode_pk = request.GET.get("barcodePk", None)
    queryset = PafSummaryCov.objects.filter(
        job_master__flowcell_id=pk, chromosome_pk=chromosome_pk, barcode_id=barcode_pk
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
            "chromosome_name",
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
            "chromosome_name",
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
            "chromosome_name",
            "coverage",
            "average_read_length",
            "barcode_name",
            "read_type__name",
        ]
    ].T.apply(
        lambda col: {"name": col[3] + " " + col[4] + " " + col[0], "data": col[2]}
    )
    categories = np.unique(df["chromosome_name"].values)
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
            referenceName=F("reference_name"),
            referenceId=F("reference_pk"),
            barcodeName=F("barcode_name"),
            barcodeId=F("barcode_id"),
            chromosomeName=F("chromosome_name"),
            chromosomeId=F("chromosome_pk"),
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
    Get the per genome (i.e chromosome level) coverage (for the top column charts)
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
    query = PafSummaryCov.objects.filter(job_master__flowcell_id=flowcell_pk).exclude(
        job_master__job_type_id=16
    )
    categories = query.values_list("chromosome_name", flat=True).distinct()
    queryset = list(
        query
        .values("chromosome_name", "barcode_name", "reference_line_length")
        .annotate(
            Sum("total_yield"),
            Avg("average_read_length"),
        )
    )
    a = {f"{q['chromosome_name']}--{q['barcode_name']}": q["total_yield__sum"]/q["reference_line_length"] for q in queryset}
    cov_series = [{"name": barcode_name, "data": [a.get(f"{cat}--{barcode_name}", 0) for cat in categories]} for
              barcode_name in query.values_list("barcode_name", flat=True).distinct()]
    a = {f"{q['chromosome_name']}--{q['barcode_name']}": q["average_read_length__avg"] for q in queryset}
    read_len_series = [{"name": barcode_name, "data": [a.get(f"{cat}--{barcode_name}", 0) for cat in categories]} for
                  barcode_name in query.values_list("barcode_name", flat=True).distinct()]
    results = {
        "coverageData": cov_series,
        "avgRLData": read_len_series,
        "categories": categories
    }
    return Response(results, status=status.HTTP_200_OK)
