import json

import numpy as np
import pandas as pd
from django.db.models import F, Sum, Avg
from natsort import natsorted
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from alignment.models import PafSummaryCov
from alignment.utils import _human_key, get_alignment_result_dir, get_or_create_array
from minknow_data.models import Flowcell
from reads.models import Barcode, JobMaster
from reference.models import ReferenceLine


def sampler(array, array_length=5000):
    """
    If the array is too large > 50000 points, take a sample using steps until we have 50000 or less points
    :param array: A list of tuple values for expected benefit steps
    :type array: list
    :return: a smaller sampled numpy array
    """
    # Get the step size for slicing, getting the array to 50000 points
    step_size = int(len(array) / array_length)
    # Slice is step size
    sampled_array = array[0::step_size]
    return sampled_array


@api_view(["GET"])
def get_cnv_barcodes(request, pk):
    """
    Get the CNV barcodes
    Parameters
    ----------
    request: rest_framework.request.Request
        The rest framework request object
    pk: int
        The primary key of the flowcell

    Returns
    -------
    None
    """
    flowcell = Flowcell.objects.get(pk=pk, owner_id=request.user.id)
    runs = flowcell.runs.all().values_list("id", flat=True)
    jobs = JobMaster.objects.filter(run_id__in=runs)
    barcode_results = set()
    for job in jobs:
        barcodes = set(
            PafSummaryCov.objects.filter(job_master=job).values_list(
                "barcode_name", "barcode_id"
            )
        )
        barcode_results.update(barcodes)

    return Response(
        natsorted(barcode_results, key=lambda x: x[0]), status=status.HTTP_200_OK
    )


@api_view(["GET"])
def coverage_detail(
    request, task_id, barcode_id, read_type_id, chromosome_id, start, end
):
    """
    Fetch data for coverage detail charts
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
    start: int
        Position from which you want data on a chromosome
    end: int
        Position from which you want data to end on a chromosome.

    Returns
    -------

    """

    # Get the barcode name so we can recover all the results for this run.
    barcode = Barcode.objects.get(pk=barcode_id)
    job_master = JobMaster.objects.get(pk=task_id)
    # todo i guess this should be different (add a run to the drop downs)
    run = job_master.flowcell.runs.last()
    contig = ReferenceLine.objects.get(pk=chromosome_id)
    folder_dir = get_alignment_result_dir(run.runid, create=False)
    array_path = get_or_create_array(
        folder_dir, contig_name=contig.line_name, barcode_name=barcode.name
    )
    if not array_path.exists():
        return Response(status=status.HTTP_404_NOT_FOUND)
    mem_map = np.load(array_path, mmap_mode="r")
    # TODO limits to just one reference here when fetching length, could be an issue
    #  in the future displaying multiple references on a plot
    length = contig.chromosome_length
    x_coords = np.arange(int(start), int(end) + 10, 10, dtype=int)
    if not x_coords[0]:
        x_coords[0] += 1

    start, end = int(start) // 10, np.ceil(int(end) / 10).astype(int)
    factor = int((end - start) / 1800)
    factor = 1 if factor < 1 else factor
    x_coords = x_coords[::factor]
    rejected_new_data = mem_map[1, start:end:factor]
    sequenced_new_data = mem_map[0, start:end:factor]
    seq_coords = list(zip(x_coords, sequenced_new_data))
    rej_coords = list(zip(x_coords, rejected_new_data))
    # sum_to_check = mem_map.sum()
    new_data = {"sequenced": seq_coords, "unblocked": rej_coords}

    return Response(
        {"newChartData": new_data, "refLength": length, "sumToCheck": 0},
        status=status.HTTP_200_OK,
    )


@api_view(["GET"])
def coverage_master(request, task_id, barcode_id, read_type_id, chromosome_id):
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
    # Get the barcode name so we can recover all the results for this run.
    pixels = 3600
    barcode = Barcode.objects.get(pk=barcode_id)
    job_master = JobMaster.objects.get(pk=task_id)
    # todo i guess this should be different (add a run to the drop downs)
    run = job_master.flowcell.runs.last()
    contig = ReferenceLine.objects.get(pk=chromosome_id)
    folder_dir = get_alignment_result_dir(run.runid, create=False)
    array_path = get_or_create_array(
        folder_dir, contig_name=contig.line_name, barcode_name=barcode.name
    )
    if not array_path.exists():
        return Response(status=status.HTTP_404_NOT_FOUND)
    mem_map = np.load(array_path, mmap_mode="r")
    # TODO limits to just one reference here when fetching length, could be an issue
    #  in the future displaying multiple references on a plot
    length = contig.chromosome_length
    step = int((length / 10) / pixels)
    step = 1 if step < 1 else step
    data = mem_map[::step]
    rejected_new_data = mem_map[1, ::step]
    sequenced_new_data = mem_map[0, ::step]
    x_coords = np.arange(0, length + 10, 10)[::step]
    seq_coords = list(zip(x_coords, sequenced_new_data))
    rej_coords = list(zip(x_coords, rejected_new_data))
    new_data = {"sequenced": seq_coords, "unblocked": rej_coords}
    sum_to_check = mem_map.sum()
    print(sum_to_check)
    return Response(
        {"newChartData": new_data, "refLength": length, "sumToCheck": sum_to_check},
        status=status.HTTP_200_OK,
    )


@api_view(["GET"])
def cnv_chart(request, pk: int, barcode_pk: int, expected_ploidy: int):
    """
    get the xy scatter coordinates of the CNV ploidy chart
    Parameters
    ----------
    request: rest_framework.Request.request
        Django rest framework object
    pk: int
        The primary key of the flowcell
    barcode_pk: int
        Primary key of the barcode to look at
    expected_ploidy: int
        The expected ploidy, should be 2

    Returns
    -------
    dict of lists
        Sum to check, to redraw
        new chart data in list of tuple (xy coordinates)
    """
    expected_ploidy = int(expected_ploidy)
    flowcell = Flowcell.objects.get(pk=pk)
    barcode = Barcode.objects.get(pk=barcode_pk).name
    runs = flowcell.runs.all().values_list("id", flat=True)
    jobs = JobMaster.objects.filter(run_id__in=runs)
    array_path_me_baby = {}
    result_me_baby = {}
    # reads_per_bin = 10
    for job in jobs:
        genome_length = job.reference.length
        result_dir = get_alignment_result_dir(job.run.runid, create=False)
        # bin size
        cumsum_chromosome_lengths = np.cumsum(
            [0]
            + list(
                map(
                    int,
                    job.reference.reference_lines.values_list(
                        "chromosome_length", flat=True
                    ),
                )
            )
        )
        # we need this below in order to move the contig coordinates the correct amount along the graph
        chromo_name_to_length = dict(
            zip(
                job.reference.reference_lines.values_list("line_name", flat=True),
                cumsum_chromosome_lengths,
            )
        )
        total_map_starts = 0
        array_path_me_baby[barcode] = natsorted(
            list(result_dir.rglob(f"{barcode}/*/*cnv_bins.npy")),
            key=lambda x: x.parts[-2],
        )
        # get the total mapping starts per barcode
        for contig_array_path in array_path_me_baby[barcode]:
            if contig_array_path.parts[-2] == "chrM":
                continue
            contig_array_mmap = np.load(contig_array_path, mmap_mode="r")
            total_map_starts += contig_array_mmap[0].sum()
        bin_size = 50000
        for contig_array_path in array_path_me_baby[barcode]:
            contig_name = contig_array_path.parts[-2]
            if contig_name == "chrM":
                continue
            bin_slice = np.ceil(bin_size / 10).astype(int)
            contig_array = np.load(contig_array_path)
            new_bin_values = np.fromiter(
                (
                    contig_array[0][start : start + bin_slice].sum()
                    for start in range(0, contig_array[0].size + 1, bin_slice)
                ),
                dtype=np.float64,
            )
            median_bin_value = np.median(new_bin_values)
            x_coords = range(
                0 + chromo_name_to_length[contig_name],
                (contig_array[0].size + 1) * 10
                + chromo_name_to_length[contig_name],
                bin_slice * 10,
            )
            binned_ploidys = np.nan_to_num(
                new_bin_values / median_bin_value * expected_ploidy,
                nan=0,
                posinf=0,
            ).round(decimals=5)
            
            result_me_baby[contig_name] = np.array(
                list(zip(
                    x_coords, binned_ploidys
                ))
            )[binned_ploidys!=0]
            # algo_c = rpt.KernelCPD(kernel="linear", min_size=10).fit(
            #     binned_ploidys
            # )
            # # guesstimated parameter for how big a shift needs to be.
            # penalty_value = 3
            #
            # # The value of my_bkps is the coordinate in the sample data where a break occurs.
            # # This therefore needs to be converted back to the BIN value where this occured.
            # # thus for each point in the my_bkps we need to grab the row from the data table.
            #
            # my_bkps = algo_c.predict(pen=penalty_value)
        result_me_baby["plotting_data"] = chromo_name_to_length
        return Response(result_me_baby)


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
    df = pd.DataFrame.from_records(queryset)
    df2 = (
        df.groupby(
            by=[
                "barcode_name",
                "chromosome_name",
                "job_master_id",
                "reference_line_length",
                "reference_name",
            ]
        )
        .sum()
        .reset_index()
    )
    # df2['average_read_length'] = (df2['total_yield'] / df2['read_count']).astype(int)
    # df2['coverage'] = (df2['total_yield'] / df2['reference_line_length']).astype(float)

    result = {
        # "data": queryset,
        "data": df2.to_dict(orient="records"),
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

    # ToDo: This is returning multiple barcodes when really we want to return just one.
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
    return Response({"references": references}, status=status.HTTP_200_OK)


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
    categories = list(query.values_list("chromosome_name", flat=True).distinct())
    categories.sort(key=_human_key)
    queryset = list(
        query.values(
            "chromosome_name", "barcode_name", "reference_line_length"
        ).annotate(
            Sum("total_yield"), Avg("average_read_length"),
        )
    )
    a = {
        f"{q['chromosome_name']}--{q['barcode_name']}": q["total_yield__sum"]
        / q["reference_line_length"]
        for q in queryset
    }
    cov_series = [
        {
            "name": barcode_name,
            "data": [a.get(f"{cat}--{barcode_name}", 0) for cat in categories],
        }
        for barcode_name in query.values_list("barcode_name", flat=True).distinct()
    ]
    a = {
        f"{q['chromosome_name']}--{q['barcode_name']}": q["average_read_length__avg"]
        for q in queryset
    }
    read_len_series = [
        {
            "name": barcode_name,
            "data": [a.get(f"{cat}--{barcode_name}", 0) for cat in categories],
        }
        for barcode_name in query.values_list("barcode_name", flat=True).distinct()
    ]
    results = {
        "coverageData": cov_series,
        "avgRLData": read_len_series,
        "categories": categories,
    }
    return Response(results, status=status.HTTP_200_OK)
