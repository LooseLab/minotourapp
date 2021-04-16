"""
views.py
"""
import math
import tarfile

import numpy as np
import pandas as pd
from django.contrib.auth.decorators import login_required
from django.db.models import Sum, Q
from django.http import HttpResponse
from django.shortcuts import render
from django.utils import timezone
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.decorators import api_view
from rest_framework.response import Response
from weasyprint import HTML, CSS

from alignment.models import PafSummaryCov
from metagenomics.models import (
    CentrifugeOutput,
    MappingResult,
    Metadata,
    SankeyLink,
    DonutData,
    MappingTarget,
    EstimatedAbundance, UncertaintyProbability,
)
from metagenomics.utils import calculate_proportion_for_table, get_metagenomics_data
from minknow_data.models import Flowcell
from reads.models import JobMaster, JobType

pd.options.mode.chained_assignment = None


@api_view(["GET", "POST"])
def export_report(request, pk):
    """
    Export the reach back report to the user upon request
    Parameters
    ----------
    request: rest_framework.request.Request
        The request object
    pk: int
        Flowcell primary key

    Returns
    -------
    A PDF blob to be downloaded
    """
    flowcell = Flowcell.objects.get(pk=pk)
    task = JobMaster.objects.filter(flowcell_id=pk, job_type_id=10).last()
    run = flowcell.runs.first()
    metagenomics_run = True if task else False
    alignment_data = PafSummaryCov.objects.filter(job_master__flowcell_id=pk).exclude(job_master__job_type__id=16).values()
    report_info = {
        "run_id": run.name,
        "read_count": run.summary.read_count,
        "yield": run.summary.total_read_length,
        "avg_read_length": run.summary.avg_read_length,
        "run_date": run.summary.first_read_start_time.strftime("%d/%m/%Y"),
        "run_duration": str(
            run.summary.last_read_start_time - run.summary.first_read_start_time
        ),
        "alignment_data": bool(alignment_data),
        "metagenomics_data": metagenomics_run
    }
    if metagenomics_run:
        metadata = Metadata.objects.get(task=task)
        meta_dict = {
            "unclassified": metadata.unclassified,
            "classified": metadata.classified,
            "percent_classified": round(
                metadata.classified / run.summary.read_count * 100, 2
            ),
            "percent_unclassified": round(
                metadata.unclassified / run.summary.read_count * 100, 2
            ),
        }
        report_info.update(meta_dict)
    if alignment_data:
        alignment_dict = {
            "alignment_values": alignment_data
        }
        report_info.update(alignment_dict)
    print(request.data["basecalledSummaryHTML"])
    print(request.data["liveEventSummaryHTML"])
    report_info.update(request.data)
    cent_df = get_metagenomics_data(task, flowcell)
    cent_df.to_csv(f"/tmp/meta_csv_{run.name}.csv.gz", encoding="utf-8")
    HTML(
        string=render(
            request, "metagenomics/report.html", context={"report_data": report_info}
        ).getvalue()
    ).write_pdf(
        f"/tmp/{task.id}_run_{run.name}_report.pdf",
        stylesheets=[
            CSS("web/static/web/css/report.css"),
            CSS("web/static/web/libraries/bootstrap-4.5.0-dist/css/bootstrap.css"),
        ],
    )
    tar_file_path = f"/tmp/{task.id}_run_{run.name}_report.tar.gz"
    tar_file_name = f"{task.id}_run_{run.name}_report.tar.gz"
    with tarfile.open(tar_file_path, "w:gz") as tar:
        try:
            tar.add(f"/tmp/{task.id}_run_{run.name}_report.pdf", recursive=False)
            tar.add(f"/tmp/meta_csv_{run.name}.csv.gz", recursive=False)
        except FileNotFoundError as e:
            print("file not found")

    with open(tar_file_path, "rb") as fh:
        response = HttpResponse(fh.read(), content_type="application/zip")
        response[
            "Content-Disposition"
        ] = f"attachment; filename={tar_file_name}"
        response["x-file-name"] = tar_file_name

    return response


@api_view(["GET"])
def get_abundance_table_html(request, pk):
    """

    Parameters
    ----------
    request: rest_framework.request.Request
        The AJAX request
    pk: int
        The primary key of the flowcell that the metagenomics task is on

    Returns
    -------
    str
        Pre formated HTML string
    """
    task = JobMaster.objects.filter(job_type_id=10, flowcell_id=pk).last()
    targets = MappingTarget.objects.filter(target_set=task.target_set)
    uncertainties = UncertaintyProbability.objects.filter(task=task).values_list("lower_ci_value", "upper_ci_value")
    df = pd.DataFrame.from_records(
        EstimatedAbundance.objects.filter(task=task).values(
            "tax_id", "abundance", "name"
        )
    )
    df["Lower Conf. bound"] = 0
    df["Higher Conf. bound"] = 0
    for ind, (low, high) in enumerate(uncertainties):
        print(low)
        df.loc[ind, "Lower Conf. bound"] = low
        df.loc[ind, "Higher Conf. bound"] = high
    df = df.set_index("tax_id")
    for target in targets:
        if int(target.tax_id) not in df.index:
            df.loc[int(target.tax_id)] = [0, target.species, 0, 0]

    df = df.reset_index()
    return Response(
        df.to_html(classes="table table-striped", border=0, justify="left")
    )


def check_run_is_legit(flowcell, task):
    """
    Check run has enough data, and other metrics are met

    Parameters
    ----------
    flowcell: minknow_data.models.Flowcell
        Flowcell ORM
    task: reads.models.JobMaster
        JobMaster ORM object

    Returns
    -------
    (bool, list of str)
        bool if the run is legit or not, if not the reason it isn't. Default is a blank string.
    """
    queryset = Metadata.objects.filter(task=task).last()
    read_threshold = 50000
    classified_ratio = 20
    num_reads_check = task.read_count <= read_threshold
    # data_generated_hour_rate_check =
    twenty_classified_check = (
        queryset.classified / task.read_count * 100
    ) <= classified_ratio
    reason = np.array(
        [
            f"Number reads is below threshold {read_threshold}",
            f"Ratio of classified to unclassified is below {classified_ratio}",
        ]
    )
    success = [num_reads_check, twenty_classified_check]
    reason = reason[success]
    return any(success), reason.tolist()


@api_view(["GET"])
def super_simple_alert_list(request, pk, barcode_name):
    """
    SUper simple alert table. Checks if run meets conditions for trustworthyness and if the alerts
    Parameters
    ----------
    request: rest_framework.request.Request
        Ajax request body

    Returns
    -------

    """
    flowcell = Flowcell.objects.get(pk=pk)
    task = (
        JobMaster.objects.filter(flowcell__id=pk, job_type__name="Metagenomics")
        .order_by("id")
        .last()
    )
    task_id = int(task.id)
    queryset = MappingResult.objects.filter(
        task__id=task_id, barcode_name=barcode_name.replace("_", " ")
    ).values()
    results_df = pd.DataFrame.from_records(queryset)
    results_df["detected"] = np.where(results_df["num_mapped"] > 0, True, False)
    results_df["percent"] = (
        results_df["num_mapped"].div(task.read_count).mul(100).round(2).astype(str)
    )
    results_df["reason"] = (
        "Species "
        + results_df["species"]
        + " identified at "
        + results_df["percent"]
        + "% of reads sequenced."
    )
    run_status_bool, run_reasons = check_run_is_legit(flowcell, task)
    alert = results_df["detected"].any()
    if not alert and not run_status_bool:
        run_reasons = ["No target species detected and Run is of sufficient quality."]
    return Response(
        [
            {
                "Alert": alert,
                "run_status": run_status_bool,
                "run_status_reasons": run_reasons,
                "alert_reasons": results_df[results_df["detected"]][
                    "reason"
                ].values.tolist(),
            }
        ],
        status=status.HTTP_200_OK,
    )


def alert_level(col):
    """
    Returns the highest alert level for a barcode
    :param col: The dataframe column that we are applying the function to
    :return:
    """
    alert_levels = {"num_matches": 1, "num_mapped": 2, "red_reads": 3}
    series = np.where(col > 0, alert_levels[col.name], 0)
    return series


def rounddown(x, sigfigs=1):
    """
    Round down to one significant figure
    :param x: The number to round
    :param sigfigs: The number of significant figures to round down to
    :return:
    """
    # towards -inf
    exponent = math.floor(
        math.log10(math.copysign(x, 1))
    )  # we don't want to accidentally try and get an
    #  imaginary log (it won't work anyway)
    mantissa = x / 10 ** exponent  # get full precision mantissa
    # change floor here to ceil or round to round up or to zero
    mantissa = math.floor(mantissa * 10 ** (sigfigs - 1)) / 10 ** (
        sigfigs - 1
    )  # round mantissa to sigfigs
    return int(mantissa * 10 ** exponent)


@api_view(["GET"])
def centrifuge_metadata(request):
    """
    :purpose: return the metadata displayed at the top of the page relating to the task data being currently visualised

    :author: Rory

    :param request: A django rest framework request body object, contains the flowcell id as a query parameter
    :return: A Response object, containing a list with the four data objects, Reads sequenced, Reads Classified,
    Classified and the Runtime
    """
    # The flowcell id for the flowcell that the fastq data came from, default to False if not present
    flowcell_id = request.GET.get("flowcellId", False)
    # The ids of the JobMaster entries for this flowcell
    job_master = (
        JobMaster.objects.filter(
            flowcell__id=flowcell_id, job_type__name="Metagenomics"
        )
        .order_by("id")
        .last()
    )

    if not job_master:
        return Response(
            "No Centrifuge tasks. This is not the API you are looking for....",
            status=404,
        )

    # If there is no MetaGenomicsMeta object return an error
    try:
        queryset = Metadata.objects.filter(task__id=job_master.id).last()
    except Metadata.DoesNotExist as e:
        return Response(e, status=404)

    if not queryset:
        return Response(
            "No Metadata. This is not the task you are looking for....", status=404
        )

    print(job_master.flowcell.number_reads)
    # Number of classified reads
    reads_actually_classified = queryset.classified
    # Number of unclassified reads
    reads_not_classified = queryset.unclassified
    # get the number of reads
    number_of_reads = Flowcell.objects.get(pk=flowcell_id).number_reads
    # Get the read count
    reads_class = job_master.read_count
    # Percentage of reads classified
    percentage = round(reads_class / number_of_reads * 100, 2)
    # Classified as a percentage
    class_percent = round(reads_actually_classified / reads_class * 100, 2)
    # Unclassified as a percentage
    unclass_percent = round(reads_not_classified / reads_class * 100, 2)
    # Get the start time of the task, removing the timezone info
    start_time = queryset.start_time.replace(tzinfo=None)
    # Get the current time, removing the timezone info
    current_time = timezone.now().replace(tzinfo=None)
    # If the request is not for a finished run, subtract the start time from the current
    # time to give the time taken so far
    if not queryset.finish_time:
        runtime = current_time - start_time
        runtime = str(runtime)
    # If the run is finished, use the final time taken saved into the database
    else:
        runtime = queryset.finish_time

    # this bool tells the javascript whether not to place the metadata in the validation
    #  set panel or the analysis results panel
    print(job_master.target_set)
    if job_master.target_set == "-1" or job_master.target_set == None:
        validation_table_present = False
    else:
        print("Showing validation table")
        validation_table_present = True

    return_list = {
        "validation": validation_table_present,
        "result": [
            {"key": "Reads Sequenced: ", "value": number_of_reads},
            {
                "key": "Reads Analysed: ",
                "value": str(reads_class) + " (" + str(percentage) + "%)",
            },
            {
                "key": "Reads Classified: ",
                "value": str(reads_actually_classified)
                + " ("
                + str(class_percent)
                + "%)",
            },
            {
                "key": "Reads Unclassified: ",
                "value": str(reads_not_classified) + " (" + str(unclass_percent) + "%)",
            },
            {"key": "Runtime: ", "value": runtime},
        ],
    }

    return Response(return_list)


@api_view(["GET"])
def centrifuge_sankey(request):
    """
    :purpose: Query the database for the sankeyLink data, return the top 50 Lineages
    :author: Rory
    :param request: (obj) Django rest framework object, with query params, speciesLimit -
    The number of species to Return
    flowcell_id, the Flowcell id for the flowcell that provided the reads
    :return: (obj) A Response obj containing a dict with data for the sankey diagram
    """

    # Get the number of species to visualise , defaulting to 50 if not present
    species_limit = request.GET.get("speciesLimit", 20)
    # Get the flowcell ID , defaulting to False if not present
    flowcell_id = request.GET.get("flowcellId", False)
    # Selected barcode, default all reads
    selected_barcode = request.GET.get("barcode", "All reads")
    # The most up to date task_id
    task_id = (
        JobMaster.objects.filter(
            flowcell__id=flowcell_id, job_type__name="Metagenomics"
        )
        .order_by("id")
        .last()
        .id
    )

    # ## Get the jobMaster for calculate sankey
    been_run = (
        JobMaster.objects.filter(
            flowcell__id=flowcell_id, job_type__name="CalculateMetagenomicsSankey"
        )
        .order_by("id")
        .last()
    )
    # If the sankey task has been run
    if been_run:
        has_been_run = True
    else:
        return Response({"message": "No Sankey data has been found", "run": False})
    # Get the sankey links objects for this run
    queryset = SankeyLink.objects.filter(
        task__id=task_id, barcode_name=selected_barcode
    ).values()
    # If the queryset is empty, return an empty object
    if not queryset:
        return Response({"message": "No Sankey data has been found", "run": False})
    # Create the dataframe from the results of the database query
    source_target_df = pd.DataFrame(list(queryset))
    # If the database is empty return an empty list
    if source_target_df.empty:
        return Response({"message": "No Sankey data has been found", "run": False})
    # get a subset df of all the species rows
    temp_species_df = source_target_df[
        source_target_df["target_tax_level"] == "species"
    ]
    # get species limit (default 50) of the largest species
    temp_species_df = temp_species_df.nlargest(species_limit, ["value"])
    # get the values for those species
    source_target_df = source_target_df[
        source_target_df["tax_id"].isin(temp_species_df["tax_id"])
    ]
    source_target_df.set_index(["source", "target"], inplace=True)
    # Group by rows where the source and target are the same
    st_gb = source_target_df.groupby(["source", "target"])
    # Replace the value columns on dataframe with the sum of all the values of identical source target rows
    source_target_df["value"] = st_gb["value"].sum()

    source_target_df.reset_index(inplace=True)
    # Drop all duplicate rows, only need one new entry
    source_target_df.drop_duplicates(["source", "target"], keep="first", inplace=True)
    # Drop any rows where the source and the target are the same and don't keep them
    source_target_df = pd.concat(
        [
            source_target_df,
            source_target_df[source_target_df["source"] == source_target_df["target"]],
        ]
    ).drop_duplicates(["source", "target"], keep=False)
    source_target_df.dropna(inplace=True)

    # Drop unnecessary columns from DataFrame
    source_target_df.drop(
        columns=["tax_id", "target_tax_level", "task_id"], inplace=True
    )
    # Set MultiIndex to group source to target
    source_target_df = pd.concat(
        [
            source_target_df,
            source_target_df[source_target_df["source"] == source_target_df["target"]],
        ]
    ).drop_duplicates(["source", "target"], keep=False)
    # source_target_df.sort_values(["value"], ascending=False, inplace=True)
    # Create the links list of dicts to return
    links = source_target_df.to_dict(orient="records")

    # ## Get the nodes ###

    # Get all the nodes values from superkingdom (ex. Bacteria) to species ("E. Coli")
    # Create a series of all possible nodes
    nodes = source_target_df["source"].append(source_target_df["target"])

    # Remove duplicates
    nodes = pd.DataFrame({"name": nodes.unique()})
    # Merge the nodes and links dataframes together
    merge = pd.merge(
        nodes, source_target_df, how="outer", right_on="source", left_on="name"
    )
    # Get just the path and name for the nodes data
    merge = merge[["name", "path"]]
    # Do the same for the targets
    merge = pd.merge(
        merge, source_target_df, how="outer", right_on="target", left_on="name"
    )
    # fill the missing values from not having the sources on the links with the values from merge
    merge["path_x"].fillna(merge["path_y"], inplace=True)
    # get the complete set of names and paths from the dataframe
    merge = merge[["name", "path_x"]]
    # Drop duplicates on name, as we need unique node names
    merge.drop_duplicates(subset="name", inplace=True, keep="first")

    merge.rename(columns={"path_x": "path"}, inplace=True)

    # Put the data in the right format
    nodes = merge.to_dict(orient="records")

    # ## Return the array # ##
    nodes = {"links": links, "nodes": nodes}

    return_dict = {"sankey": nodes, "run": has_been_run}

    return Response(return_dict, status=200)


@api_view(["GET"])
def donut_data(request):
    """
    :purpose: Create and return the data for the bottom of the page total results table or the
     data for the donut chart
    :author: Rory
    :param request: (obj) - A django rest framework request object with query params - FlowcellID: the flowcell ID
    visType  - Whether this is a a request for donut chart or results table data
    :return:
    """
    flowcell_id = request.GET.get("flowcellId", 0)
    barcode = request.GET.get("barcode", "All reads")
    # Get the most recent metagenomics job
    task = (
        JobMaster.objects.filter(
            flowcell__id=flowcell_id, job_type__name="Metagenomics"
        )
        .order_by("id")
        .last()
    )
    # Initialise a dataframe
    results_df = pd.DataFrame()
    # The tax rank filter is the taxonomic ranks in order
    tax_rank_filter = [
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
    # Get the data for each tax rank, append it to the all results results_df
    for rank in tax_rank_filter:
        data_df = pd.DataFrame(
            list(
                DonutData.objects.filter(task=task, barcode_name=barcode, tax_rank=rank)
                .values("name", "tax_rank")
                .distinct()
                .annotate(sum_value=Sum("num_matches"))
                .order_by("-sum_value")[0:10]
            )
        )
        results_df = results_df.append(data_df)
    if results_df.empty:
        return Response("No results", status=204)
    # We only need this data
    results_df = results_df[["tax_rank", "name", "sum_value"]]
    # Rename the columns again
    results_df = results_df.rename(columns={"name": "label", "sum_value": "value"})
    # Error trap if there are NaNs
    results_df.fillna("Unclassified", inplace=True)
    # Group the results into by their taxonomic rank
    gb = results_df.groupby("tax_rank")
    # A return dictionary
    return_dict = {}
    # Run through the groupby by group, adding the results to the dictionary
    for name, group in gb:
        return_dict[name] = group.to_dict(orient="records")
    # Return the results
    return Response(return_dict, status=200)


@api_view(["GET"])
def get_target_mapping(request):
    """
    Get the target species data for creating the more complex table
    :param request:
    :return:
    """
    # The id of the flowcell that produced these reads
    flowcell_id = request.GET.get("flowcellId", 0)
    # The barcode that is currently selected to be viewed on the page
    barcode = request.GET.get("barcode", "All reads")
    if flowcell_id == 0:
        return Response("Flowcell id has failed to be delivered", status=404)
    # Get the most recent jobmaster id, although there should only be one
    task = (
        JobMaster.objects.filter(
            flowcell__id=flowcell_id, job_type__name="Metagenomics"
        )
        .order_by("id")
        .last()
    )
    # If the barcode is All reads, there is always four
    queryset = MappingResult.objects.filter(task=task, barcode_name=barcode).values()
    results_df = pd.DataFrame(list(queryset))
    if results_df.empty:
        return Response(
            "No data has yet been produced for the target mappings", status=204
        )
    results_df[
        [
            "proportion_of_classified",
            "mapped_proportion_of_classified",
            "red_reads_proportion_of_classified",
        ]
    ] = pd.DataFrame(
        np.array(
            results_df.apply(
                calculate_proportion_for_table, args=(task,), axis=1
            ).values.tolist()
        )
    )
    results_df.rename(
        columns={
            "num_mapped": "Num. mapped",
            "red_reads": "Target reads",
            "red_sum_unique": "Unique Target reads",
            "species": "Species",
            "tax_id": "Tax id",
            "num_matches": "Num. matches",
            "sum_unique": "Sum. Unique",
            "mapped_proportion_of_classified": "Mapped prop. matches (%)",
            "red_reads_proportion_of_classified": "Target prop. mapped (%)",
            "proportion_of_classified": "Prop. classified (%)",
        },
        inplace=True,
    )
    ## calculate proportion mapped
    results = results_df.to_dict(orient="records")
    return_dict = {"table": results}
    return Response(return_dict)


@api_view(["GET"])
def metagenomic_barcodes(request, pk):
    """
    :param request:
    :param pk:
    :return:
    """
    flowcell_list = Flowcell.objects.filter(owner=request.user).filter(id=pk)
    metagenomics_barcodes = []
    if flowcell_list.count() > 0:
        flowcell = flowcell_list[0]
        task = (
            JobMaster.objects.filter(flowcell=flowcell, job_type__name="Metagenomics")
            .order_by("id")
            .last()
        )
        if task:
            metagenomics_barcodes = (
                CentrifugeOutput.objects.filter(task__id=task.id)
                .exclude(barcode_name="No")
                .values_list("barcode_name", flat=True)
                .distinct()
            )
    else:
        return Response("No flowcells found for this owner", status=204)
    # #### Get the barcode alert levels #### #
    task_id = (
        JobMaster.objects.filter(flowcell=flowcell, job_type__name="Metagenomics")
        .order_by("id")
        .last()
        .id
    )
    barcode_df = pd.DataFrame(
        list(MappingResult.objects.filter(task__id=task_id).values())
    )
    if not barcode_df.empty:
        calc_df = barcode_df[["num_matches", "num_mapped", "red_reads"]]
        hodf = calc_df.apply(alert_level)
        barcode_df["alert_level"] = hodf.max(axis=1)
        alert_level_results = barcode_df.groupby("barcode_name")["alert_level"].max()
        alert_level_results = alert_level_results.to_dict()
    else:
        alert_level_results = {}
    # Metagenomics barcodes is a list of the barcode_names found in this flowcell, tabs is the alert level 0-3,
    #  3 being target found
    return Response({"data": metagenomics_barcodes, "tabs": alert_level_results})


@api_view(["GET"])
def all_results_table(request, pk):
    """
    Returns data for metagenomics results table, containing all species matched
    Parameters
    ----------
    request: rest_framework.request.Request
        Request object body
    pk: int
        Primary key of the flowcell record in database

    Returns
    -------

    """
    print("STARTED")
    query_columns = [
        "barcode_name",
        "superkingdom",
        "phylum",
        "classy",
        "order",
        "family",
        "genus",
        "species",
        "num_matches",
        "proportion_of_classified",
    ]
    flowcell_id = pk
    draw = int(request.GET.get("draw", 0))
    start = int(request.GET.get("start", 0))
    length = int(request.GET.get("length", 10))
    end = start + length
    # Which column
    search_value = request.GET.get("search[value]", "")
    order_column = request.GET.get("order[0][column]", "")
    # ascending descending
    order_dir = request.GET.get("order[0][dir]", "")
    meta_task = JobMaster.objects.get(
        flowcell__id=flowcell_id, job_type__name="Metagenomics"
    )
    latest = meta_task.iteration_count
    if not search_value == "":
        cent_out_temp = (
            CentrifugeOutput.objects.filter(
                task__flowcell_id=flowcell_id, latest=latest
            )
            .filter(task__flowcell__owner=request.user)
            .filter(species__icontains=search_value)
            | CentrifugeOutput.objects.filter(
                task__flowcell_id=flowcell_id, latest=latest
            )
            .filter(task__flowcell__owner=request.user)
            .filter(genus__icontains=search_value)
            | CentrifugeOutput.objects.filter(
                task__flowcell_id=flowcell_id, latest=latest
            )
            .filter(task__flowcell__owner=request.user)
            .filter(family__icontains=search_value)
            | CentrifugeOutput.objects.filter(
                task__flowcell_id=flowcell_id, latest=latest
            )
            .filter(task__flowcell__owner=request.user)
            .filter(order__icontains=search_value)
        )

    else:
        cent_out_temp = CentrifugeOutput.objects.filter(
            task__flowcell_id=flowcell_id, latest=latest
        ).filter(task__flowcell__owner=request.user)

    if order_column:
        if order_dir == "desc":
            cent_out_temp = cent_out_temp.order_by(
                "-{}".format(query_columns[int(order_column)])
            )
        else:
            cent_out_temp = cent_out_temp.order_by(
                "{}".format(query_columns[int(order_column)])
            )
    cents = cent_out_temp.exclude(barcode_name="No").values(
        "tax_id",
        "barcode_name",
        "name",
        "num_matches",
        "proportion_of_classified",
        "superkingdom",
        "phylum",
        "classy",
        "order",
        "family",
        "genus",
        "species",
    )
    records_total = cent_out_temp.count()
    result = {
        "draw": draw,
        "recordsTotal": records_total,
        "recordsFiltered": records_total,
        "data": list(cents[start:end]),
    }
    print(result)
    return Response(result, status=status.HTTP_200_OK)


@api_view(["GET"])
def simple_target_mappings(request):
    """
        Get the target species
        :param request:
        :return:
        """
    # The id of the flowcell that produced these reads
    flowcell_id = request.GET.get("flowcellId", 0)
    # The barcode that is currently selected to be viewed on the page
    barcode = request.GET.get("barcode", "All reads")
    if flowcell_id == 0:
        return Response("Flowcell id has failed to be delivered", status=402)
    # Get the most recent jobmaster id, although there should only be one
    task = (
        JobMaster.objects.filter(
            flowcell__id=flowcell_id, job_type__name="Metagenomics"
        )
        .order_by("id")
        .last()
    )
    number_species_identified = (
        CentrifugeOutput.objects.filter(task=task, barcode_name=barcode)
        .exclude(species__contains="Unclassified")
        .aggregate(Sum("num_matches"))["num_matches__sum"]
    )
    confidence_detection_limit = round(
        1 / (1 - pow(0.01, 1 / number_species_identified))
    )
    if confidence_detection_limit > 100000:
        sig_fig = 2
    else:
        sig_fig = 1
    confidence_detection_limit = rounddown(confidence_detection_limit, sig_fig)
    # If the barcode is All reads, there is always four
    queryset = MappingResult.objects.filter(task=task, barcode_name=barcode).values()
    if not queryset:
        return Response("Metagenomics task has no validation set", status=204)
    results_df = pd.DataFrame(list(queryset))
    results_df.drop(
        columns=[
            "id",
            "mapped_proportion_of_classified",
            "num_matches",
            "red_reads_proportion_of_classified",
            "sum_unique",
            "task_id",
            "tax_id",
        ],
        inplace=True,
    )
    results_df["Detected"] = np.where(results_df["num_mapped"] > 0, True, False)
    results_df["Not detected"] = np.where(results_df["num_mapped"] == 0, True, False)
    results_df["read_count"] = task.read_count
    results_df["conf_limit"] = confidence_detection_limit
    results_df["detected_at"] = np.floor(
        results_df["read_count"] / results_df["num_mapped"]
    )
    results_df.replace([np.inf, -np.inf], 0, inplace=True)
    results_df["detected_at"].fillna(0, inplace=True)
    results_df.rename(
        columns={
            "num_mapped": "Num. mapped",
            "red_reads": "Target reads",
            "species": "Validation species",
            "tax_id": "Tax id",
        },
        inplace=True,
    )
    results = results_df.to_dict(orient="records")
    return_dict = {"table": results, "conf_detect_limit": confidence_detection_limit}
    return Response(return_dict)


@api_view(["GET"])
def get_target_sets(request):
    """
    Return the target sets for the dropdown menu on the tasks page, for either the client or the task tab
    :param request: The django rest framework request object
    :return: A list of the names of the objects
    """
    cli = request.GET.get("cli", False)
    if cli:
        api_key = request.GET.get("api_key", False)
        if not api_key:
            return Response("The API key is needed to list the Target sets", status=403)
        else:
            user_id = Token.objects.get(key=api_key).user_id
            target_sets = (
                MappingTarget.objects.filter(Q(owner_id=user_id) | Q(private=False))
                .values_list("target_set", flat=True)
                .distinct()
            )
    else:
        user_id = request.user.id
        target_sets = (
            MappingTarget.objects.filter(Q(owner_id=user_id) | Q(private=False))
            .values("target_set", "id")
            .distinct()
        )
    return Response(target_sets)


@login_required
def metagenomics_data_download(request, pk):
    """
    Send the Metagenomics data back in CSV format, used by the download button found on the metagenomics tab
    metagenomics/templates/metagenomics/visualisation.html
    :param request:
    :return:
    """
    flowcell_id = request.GET.get("flowcellId", 0)
    flowcell = Flowcell.objects.get(pk=pk)
    job_type = JobType.objects.get(name="Metagenomics")
    metagenomics_task = JobMaster.objects.get(flowcell=flowcell, job_type=job_type)
    response = HttpResponse(content_type="text/csv")
    response[
        "Content-Disposition"
    ] = "attachment; filename={}_centrifuge_output.csv".format(flowcell.name)
    cent_df = get_metagenomics_data(metagenomics_task, flowcell)
    cent_df.to_csv(response)
    return response
