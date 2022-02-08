import json
import time
from collections import defaultdict
from itertools import groupby
from typing import Tuple, List

import numpy as np
import pandas as pd
import ruptures as rpt
from django.db.models import F, Sum, Avg
from django.shortcuts import get_object_or_404
from natsort import natsorted
from pandas import Series
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.views import APIView

from alignment.models import PafSummaryCov
from alignment.utils import _human_key, get_alignment_result_dir, get_or_create_array, open_and_load_array
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
    jobs = JobMaster.objects.filter(flowcell=flowcell, job_type_id=5)
    if not jobs:
        return Response("No CNV jobs on this flowcell", status=status.HTTP_204_NO_CONTENT)
    barcode_results = set()
    for job in jobs:
        barcodes = {*PafSummaryCov.objects.filter(job_master=job).values_list(
            "barcode_name", "barcode_id", "reference_name", "job_master_id"
        )}
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
    length = contig.chromosome_length
    x_coords = np.arange(int(start), int(end) + 10, 10, dtype=int)
    if not x_coords[0]:
        x_coords[0] += 1

    start, end = int(start) // 10, np.ceil(int(end) / 10).astype(int)
    factor = int((end - start) / 1800)
    factor = 1 if factor < 1 else factor
    x_coords = x_coords[::factor]
    folder_dir = get_alignment_result_dir(
        "", run.owner.username, job_master.flowcell.name, task_id, create=False
    )
    sum_to_check = 0
    new_data = {}
    for folder in folder_dir.iterdir():
        array_path = get_or_create_array(
            folder, contig_name=contig.line_name, barcode_name=barcode.name
        )
        count_check = 0
        while not array_path.exists():
            time.sleep(1)
            if count_check == 2:
                return Response(status=status.HTTP_404_NOT_FOUND)
            count_check += 1
        mem_map = open_and_load_array(array_path)
        # TODO limits to just one reference here when fetching length, could be an issue
        #  in the future displaying multiple references on a plot

        rejected_new_data = mem_map[1, start:end:factor]
        sequenced_new_data = mem_map[0, start:end:factor]
        if not new_data:
            new_data = {"sequenced": sequenced_new_data, "unblocked": rejected_new_data}
        else:
            new_data["sequenced"] += sequenced_new_data
            new_data["unblocked"] += rejected_new_data
        sum_to_check += mem_map.sum()
    # sum_to_check = mem_map.sum()
    seq_coords = list(zip(x_coords, new_data["sequenced"]))
    rej_coords = list(zip(x_coords, new_data["unblocked"]))
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
    length = contig.chromosome_length
    step = int((length / 10) / pixels)
    x_coords = np.arange(0, length + 10, 10)[::step]
    folder_dir = get_alignment_result_dir(
        "", run.owner.username, job_master.flowcell.name, task_id, create=False
    )
    sum_to_check = 0
    new_data = {}
    for folder in folder_dir.iterdir():
        array_path = get_or_create_array(
            folder, contig_name=contig.line_name, barcode_name=barcode.name
        )
        count_check = 0
        while not array_path.exists():
            time.sleep(1)
            if count_check == 2:
                return Response(status=status.HTTP_404_NOT_FOUND)
            count_check += 1
        mem_map = open_and_load_array(array_path)
        # TODO limits to just one reference here when fetching length, could be an issue
        #  in the future displaying multiple references on a plot
        step = 1 if step < 1 else step
        data = mem_map[::step]
        rejected_new_data = mem_map[1, ::step]
        sequenced_new_data = mem_map[0, ::step]

        if not new_data:
            new_data = {"sequenced": sequenced_new_data, "unblocked": rejected_new_data}
        else:
            new_data["sequenced"] += sequenced_new_data
            new_data["unblocked"] += rejected_new_data
        sum_to_check += mem_map.sum()
    seq_coords = list(zip(x_coords, new_data["sequenced"]))
    rej_coords = list(zip(x_coords, new_data["unblocked"]))
    new_data = {"sequenced": seq_coords, "unblocked": rej_coords}
    return Response(
        {"newChartData": new_data, "refLength": length, "sumToCheck": sum_to_check},
        status=status.HTTP_200_OK,
    )


def same_dist_elems(arr):
    if len(arr)>1:
        diff = arr[1] - arr[0]
        for x in range(1, len(arr) - 1):
            if arr[x + 1] - arr[x] != diff:
                return False
        return True
    else:
        return False

@api_view(["GET"])
def cnv_detail_chart(
    request,
    pk: int,
    barcode_pk: int,
    job_pk: int,
    contig_name: str,
    exp_ploidy: int,
    pen_value: int,
    min_diff: int,
    bin_slice: int,
    median_bin_value: int,
):
    """
    Cnv Detail chart - get the starts from one barcode
    Parameters
    ----------
    request: rest_framework.request.Request
    pk: int
        Flowcell primary key
    barcode_pk: int:
        Barcode Primary key
    job_pk: int
        The best job primary key
    contig_name: str
        Name of the contig
    exp_ploidy: int
        The expected ploidy of the genome, set in the UI
    pen_value: int
        Penalty value
    min_diff: int
        The minimum difference

    Returns
    -------
    dict of list
        Points: List[Tuple[int, int]]
        plot_band: List[int]
    """
    flowcell = Flowcell.objects.get(pk=pk)
    barcode = Barcode.objects.get(pk=barcode_pk).name
    job = get_object_or_404(JobMaster, job_type_id=5, pk=job_pk)
    expected_ploidy = int(exp_ploidy)
    results = {}
    result_dir = get_alignment_result_dir(
        "", job.flowcell.owner.username, flowcell.name, job.id, create=False
    )
    array_paths = list(result_dir.rglob(f"*/{barcode}/{contig_name}/{contig_name}_cnv_bins.npy*"))
    if not array_paths:
        return Response("No arrays found for this flowcell", status=status.HTTP_404_NOT_FOUND)
    contig_array = np.array([])
    for contig_array_path in array_paths:
        if not contig_array.size:
            contig_array = open_and_load_array(contig_array_path)
        else:
            contig_array += open_and_load_array(contig_array_path)
    # total_starts = contig_array[0].sum()
    contig_length = contig_array[0].shape[0] * 10
    # reads_per_bin = 100
    # bin_slice = np.ceil(contig_length / (total_starts / reads_per_bin)).astype(int) / 10
    bin_slice = int(bin_slice)
    #print("bin slice", bin_slice)
    new_bin_values = np.fromiter(
        (
            contig_array[0][start : start + bin_slice].sum()
            for start in range(0, contig_array[0].size + 1, bin_slice)
        ),
        dtype=np.float64,
    )
    # median_bin_value = np.median(new_bin_values[new_bin_values != 0])
    x_coords = range(0, (contig_array[0].size + 1) * 10, bin_slice * 10,)
    binned_ploidys = np.nan_to_num(
        new_bin_values / int(median_bin_value) * expected_ploidy, nan=0, posinf=0,
    ).round(decimals=5)

    if "points" not in results:
        results["points"] = np.array(list(zip(x_coords, binned_ploidys)))  # [
        # binned_ploidys != 0
        # ]
    else:
        results["points"] += np.array(list(zip(x_coords, binned_ploidys)))  # [
        # binned_ploidys != 0
        # ]

    algo_c = rpt.KernelCPD(kernel="linear", min_size=int(min_diff)).fit(
        results["points"][:, 1]
    )
    # guesstimated parameter for how big a shift needs to be.
    # The value of my_bkps is the coordinate in the sample data where a break occurs.
    # This therefore needs to be converted back to the BIN value where this occured.
    # thus for each point in the my_bkps we need to grab the row from the data table.

    # This block of code catches an error where the breakpoints are all precisely spread
    # We presume this isn't true and so recalculate the breakpoints.
    # This can end up being in an infinite loop.
    # To prevent that we set a crazy max iterations (25)
    my_bkps = None
    counter = 0
    while not my_bkps:
        counter+=1
        test1 = algo_c.predict(pen=int(pen_value))
        if len(test1) > 2 and same_dist_elems(test1):
            print("Houston we have a problem.")
        else:
            my_bkps = test1
        if counter >= 25:
            mk_bkps = None
            break

    print(my_bkps)
    if my_bkps:
        band_x_coords = []
        for x in my_bkps:
            band_x_coords.append(results["points"][x - 1][0])
        my_bkps = [0] + band_x_coords + [contig_length]
    results["plot_bands"] = my_bkps
    results["contig"] = contig_name
    return Response(results, status=status.HTTP_200_OK)


@api_view(["GET"])
def cnv_chart(request, pk: int, barcode_pk: int, job_pk: int, expected_ploidy: int):
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
    job_pk: int
        The primary key of the job_master
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
    job = JobMaster.objects.get(pk=job_pk, job_type_id=5)
    array_path_me_baby = {}
    result_me_baby = {}
    reads_per_bin = 100
    genome_length = job.reference.length
    result_dir = get_alignment_result_dir(
        "", job.flowcell.owner.username, flowcell.name, job.id, create=False
    )
    # bin size
    length_list = np.array(
        natsorted(
            job.reference.reference_lines.values_list(
                "line_name", "chromosome_length"
            ),
            key=lambda x: x[0],
        )
    )
    cumsum_chromosome_lengths = [0] + np.cumsum(
        length_list[:, 1].astype(int)
    ).tolist()

    # we need this below in order to move the contig coordinates the correct amount along the graph
    chromo_name_to_length = dict(
        zip(
            natsorted(
                job.reference.reference_lines.values_list("line_name", flat=True)
            ),
            cumsum_chromosome_lengths,
        )
    )
    total_map_starts = 0
    arrays = sorted(list(result_dir.rglob(f"*/{barcode}/*/*cnv_bins.npy*")), key=lambda x: x.parts[-2])
    g = groupby(arrays, key=lambda x: x.parts[-2])
    array_path_me_baby[barcode] = dict(natsorted({key: list(value) for key, value in g}.items(), key=lambda x: x[0]))
    # get the total mapping starts per barcode
    for contig_name, contig_array_paths in array_path_me_baby[barcode].items():
        if contig_name == "chrM":
            continue
        for contig_array_path in contig_array_paths:
            contig_array_mmap = open_and_load_array(contig_array_path)
            total_map_starts += contig_array_mmap[0].sum()
        # in bases not 10 bases
    bin_size = int(genome_length / (total_map_starts / reads_per_bin))
    bin_size = np.ceil(bin_size / 10).astype(int)
    # Get the median bin values across the whole genome in order to calculate contig ploidy against the other contigs
    temp_holder_new_bin_values = []
    for contig_name, contig_array_paths in array_path_me_baby[barcode].items():
        if contig_name == "chrM":
            continue
        contig_array = np.array([])
        for contig_array_path in contig_array_paths:
            if not contig_array.size:
                contig_array = open_and_load_array(contig_array_path)
            else:
                contig_array += open_and_load_array(contig_array_path)
        new_bin_values = np.fromiter(
            (
                contig_array[0][start: start + bin_size].sum()
                for start in range(0, contig_array[0].size + 1, bin_size)
            ),
            dtype=np.float64,
        )
        temp_holder_new_bin_values.extend(new_bin_values.tolist())
    median_bin_value = np.median(temp_holder_new_bin_values)
    for contig_name, contig_array_paths in array_path_me_baby[barcode].items():
        if contig_name == "chrM":
            continue
        contig_array = np.array([])
        for contig_array_path in contig_array_paths:
            if not contig_array.size:
                contig_array = open_and_load_array(contig_array_path)
            else:
                contig_array += open_and_load_array(contig_array_path)
        new_bin_values = np.fromiter(
            (
                contig_array[0][start : start + bin_size].sum()
                for start in range(0, contig_array[0].size + 1, bin_size)
            ),
            dtype=np.float64,
        )

        x_coords = range(
            0 + chromo_name_to_length[contig_name],
            (contig_array[0].size + 1) * 10 + chromo_name_to_length[contig_name],
            bin_size * 10,
        )
        binned_ploidys = np.nan_to_num(
            new_bin_values / median_bin_value * expected_ploidy, nan=0, posinf=0,
        ).round(decimals=2)

        result_me_baby[contig_name] = np.array(
            list(zip(x_coords, binned_ploidys))
        )  # [
        #    binned_ploidys != 0
        # ]
        points = result_me_baby[contig_name].shape[0]
        desired_points = 25000
        step = np.ceil(points / desired_points).astype(int)
        #print(f"points is {points}, step is {step}")
        step = 1 if step < 1 else step
        result_me_baby[contig_name] = result_me_baby[contig_name][::step]

    result_me_baby["plotting_data"] = chromo_name_to_length
    result_me_baby["bin_slice"] = bin_size
    result_me_baby["median_bin_value"] = median_bin_value
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


class JobMasterPk(APIView):
    """
    Get the Job MAster pk in swordfish for a task on a flowcell with minimap_cnv as the task type
    Parameters
    ----------

    """
    def get(self, request, flowcell_pk):
        jm = get_object_or_404(JobMaster, flowcell__id=flowcell_pk, job_type_id=5)
        jm = jm.__dict__
        jm.pop("_state")
        return Response(jm, status=status.HTTP_200_OK)


def _make_sv_targets(row: Series) -> list:
    """
    return readfish compatible targets
    Parameters
    ----------
    row: pd.core.series.Series
        A row from the SV pickle

    Returns
    -------
    list
        List of targets for the source and target of the SV, on both + and - strand
    """
    return [
        f"{row['source_chrom']},{int(row['source_start'])},{int(row['source_end'])},+",
        f"{row['source_chrom']},{int(row['source_start'])},{int(row['source_end'])},-",
        f"{row['target_chrom']},{int(row['target_start'])},{int(row['target_end'])},+",
        f"{row['target_chrom']},{int(row['target_start'])},{int(row['target_end'])},-",
            ]
    


@api_view(["GET"])
def get_cnv_positions(request, job_master_pk, reads_per_bin, expected_ploidy, min_diff):
    """
    `get_cnv_positions` - Get the CNV positions according to ruptures for a given job_master pk of job_type 5
    Parameters
    ----------
    request: rest_framework.request.Request
        The request object provided by django rest framework
    job_master_pk: intre
        Primary key of the minimap + CNV task
    reads_per_bin: int
        The number of goal reads in a bin
    expected_ploidy: int
        The expected ploidy of the species
    min_diff: int
        The minimum difference for ruptures
    Returns
    -------
    dict of dict of list
        Dictionary of barcodes to dictionary of chromosomes with lists of putative breakpoints
    """
    job_master = get_object_or_404(JobMaster, pk=job_master_pk, job_type_id=5)
    genome_length = int(job_master.reference.length)
    folder_dir = get_alignment_result_dir(
        "", request.user.get_username(), job_master.flowcell.name, job_master_pk, create=False
    )
    array_path_me_baby = defaultdict(dict)
    # get arrays pathlib paths
    # Matts pandas code for SV stuff
    df = pd.DataFrame()
    folder_dir = get_alignment_result_dir(
        "", request.user.get_username(), job_master.flowcell.name, job_master.id, create=False
    )
    pickle_paths = list(folder_dir.rglob("SV_data.pickle"))
    if pickle_paths:
        for pickle_path in pickle_paths:
            df = pd.concat([df, pd.read_pickle(pickle_path)])
        df = df[(df["source_coord"] != 1) & (df["target_coord"] != 1)]
        barcode_lookup = {b.id: b.name for b in Barcode.objects.filter(id__in=np.unique(df["barcode_id"].values))}
        df["barcode_name"] = df["barcode_id"].map(barcode_lookup)
        df = df.groupby("barcode_id", as_index=False).apply(
            lambda x: x[np.abs(x.source_coord - x.source_coord.mean()) <= (3 * x.source_coord.std())]).droplevel(0)
        df = df[["barcode_name", "source_chrom", "source", "target_chrom", "target"]]
        df = df[df["source_chrom"] != "chrM"]

    # back to cool CNV stuff
    arrays = natsorted(list(folder_dir.rglob(f"*/*/*cnv_bins.npy*")), key=lambda x: (x.parts[-3], x.parts[-2]))
    # group them by the barcode and contig name
    g = groupby(arrays, key=lambda x: (x.parts[-3], x.parts[-2]))
    # dict to store metadata return values
    barcode_info = defaultdict(dict)
    for (barcode, contig), group in g:
        array_path_me_baby[barcode][contig] = list(group)
    # For each barcode
    for barcode in array_path_me_baby:
        if barcode == "unclassified":
            continue
        # hold the values for each selected region of bins
        new_bin_values_holder = {}
        total_map_starts = 0
        # get the total number of mapped reads we have starts
        for contig_name, contig_array_paths in array_path_me_baby[barcode].items():
            if contig_name == "chrM":
                continue
            for contig_array_path in contig_array_paths:
                contig_array_mmap = open_and_load_array(contig_array_path)
                total_map_starts += contig_array_mmap[0].sum()
        # barcode_info[barcode]["total_mapped_starts"] = total_map_starts
        # calculate bin size in number of bins
        bin_size = int(genome_length / (total_map_starts / int(reads_per_bin)))
        bin_size = np.ceil(bin_size / 10).astype(int)
        # barcode_info[barcode]["bin_size"] = bin_size * 10 # bin size in bases
        barcode_info[barcode].update(
            {
                "name": barcode,
                "control": False,
                "min_chunks": 0,
                "max_chunks": 4,
                "single_on": "stop_receiving",
                "single_off": "unblock",
                "multi_on": "stop_receiving",
                "multi_off": "unblock",
                "no_seq": "proceed",
                "no_map": "proceed",
                "targets": []
            }
        )
        # Get the median bin values across the whole genome in order to calculate contig ploidy against the other contigs
        temp_holder_new_bin_values = []
        for contig_name, contig_array_paths in array_path_me_baby[barcode].items():
            if contig_name == "chrM":
                continue
            contig_array = np.array([])
            # this is something to do with run_ids as in combine the values across different runs?
            for contig_array_path in contig_array_paths:
                if not contig_array.size:
                    contig_array = open_and_load_array(contig_array_path)
                else:
                    contig_array += open_and_load_array(contig_array_path)
            # Create the new bin values by summing the slices of the array
            new_bin_values = np.fromiter(
                (
                    contig_array[0][start: start + bin_size].sum()
                    for start in range(0, contig_array[0].size + 1, bin_size)
                ),
                dtype=np.float64,
            )
            new_bin_values_holder[contig_name] = (new_bin_values, contig_array.size)
            temp_holder_new_bin_values.extend(new_bin_values.tolist())
        median_bin_value = np.median(temp_holder_new_bin_values)
        break_point_targets = []
        if bin_size * 10 < 1e6:
            for contig_name, (new_bin_values, num_bins) in new_bin_values_holder.items():
                binned_ploidys = np.nan_to_num(
                    new_bin_values / median_bin_value * int(expected_ploidy), nan=0, posinf=0,
                ).round(decimals=2)
                x_coords = list(range(0, (num_bins + 1) * 10, bin_size * 10, ))
                algo_c = rpt.KernelCPD(kernel="linear", min_size=int(min_diff)).fit(
                    binned_ploidys
                )
                # guesstimated parameter for how big a shift needs to be.
                # The value of my_bkps is the coordinate in the sample data where a break occurs.
                # This therefore needs to be converted back to the BIN value where this occured.
                # thus for each point in the my_bkps we need to grab the row from the data table.

                # This block of code catches an error where the breakpoints are all precisely spread
                # We presume this isn't true and so recalculate the breakpoints.
                # This can end up being in an infinite loop.
                # To prevent that we set a crazy max iterations (25)
                my_bkps = None
                counter = 0
                while not my_bkps:
                    counter += 1
                    test1 = algo_c.predict(pen=int(3))
                    if len(test1) > 2 and same_dist_elems(test1):
                        print("Houston we have a problem.")
                    else:
                        my_bkps = test1
                    if counter >= 25:
                        mk_bkps = None
                        break
                if my_bkps:
                    for x in my_bkps:
                        # Add contig name, start - 0.5 * bin width, start + 0.5 * bin_width, strand
                        # TODO make this a smart f string thing
                        break_point_targets.extend([f"{contig_name},{int(x_coords[x - 1] - 0.5 * bin_size*10)},{int(x_coords[x - 1] + 0.5 * bin_size*10)},+",f"{contig_name},{int(x_coords[x - 1] - 0.5 * bin_size*10)},{int(x_coords[x - 1] + 0.5 * bin_size*10)},-"])
            barcode_info[barcode]["targets"] = break_point_targets

            if not df.empty:
                half_bin_size = 0.5 * bin_size * 10
                # filthy but loosing will to live
                barcode_df = df[df["barcode_name"] == barcode]
                barcode_df["source_start"] = np.maximum(barcode_df["source"] - half_bin_size, 0)
                barcode_df["source_end"] = barcode_df["source"] + half_bin_size
                barcode_df["target_start"] = np.maximum(barcode_df["target"] - half_bin_size, 0)
                barcode_df["target_end"] = barcode_df["target"] + half_bin_size
                targets_list = barcode_df.apply(_make_sv_targets, axis=1)
                for targets in targets_list.values.tolist():
                    barcode_info[barcode]["targets"].extend(targets)
    return Response(barcode_info, status=status.HTTP_200_OK)


@api_view(["GET", "HEAD"])
def sv_table_list(request, flowcell_pk):
    """
    HEAD whether a given alignment run has a pickle of SVs and if it does allow clients to GET it
    Parameters
    ----------
    request: rest_framework.request.Request
        Request object
    Returns
    -------

    """

    if request.method == "HEAD":
        headers = {"x-data": False}
        job_masters = JobMaster.objects.filter(flowcell_id=flowcell_pk, job_type_id=5)
        for job_master in job_masters:
            folder_dir = get_alignment_result_dir(
                "", request.user.get_username(), job_master.flowcell.name, job_master.id, create=False
            )
            if list(folder_dir.rglob("SV_data.pickle")):
                headers = {"x-data": True}
        return Response(headers=headers)
    elif request.method == "GET":
        df = pd.DataFrame()
        job_masters = JobMaster.objects.filter(flowcell_id=flowcell_pk, job_type_id=5)
        for job_master in job_masters:
            folder_dir = get_alignment_result_dir(
                "", request.user.get_username(), job_master.flowcell.name, job_master.id, create=False
            )
            pickle_paths = list(folder_dir.rglob("SV_data.pickle"))
            for pickle_path in pickle_paths:
                df = pd.concat([df, pd.read_pickle(pickle_path)])
        df = df[(df["source_coord"] != 1) & (df["target_coord"] != 1)]
        barcode_lookup = {b.id: b.name for b in Barcode.objects.filter(id__in=np.unique(df["barcode_id"].values))}
        df["barcode_name"] = df["barcode_id"].map(barcode_lookup)
        df = df[df["source_chrom"] != "chrM"]
        df["source"] = df["source_chrom"] + " " + df["source"].astype(str)
        df["target"] = df["target_chrom"] + " " + df["target"].astype(str)
        df["read_count"] = df["source_coord"]
        df = df.drop(columns=["barcode_id", "source_chrom", "target_chrom", "source_coord", "target_coord"])
        data = df.to_dict(orient="records")
        return Response({"data": data})
