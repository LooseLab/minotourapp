"""
views.py
"""
from rest_framework.response import Response
from rest_framework.decorators import api_view
from centrifuge.models import CentrifugeOutput, MappingResult, LineageValue, Metadata, \
    SankeyLink, MappingResultsBarcoded, DonutData
from django.utils import timezone
from django.http import JsonResponse
import pandas as pd
from jobs.models import JobMaster
from reads.models import Flowcell

pd.options.mode.chained_assignment = None


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
    job_master = JobMaster.objects.filter(flowcell__id=flowcell_id, job_type__name="Metagenomics").order_by("id").last()

    if not job_master:
        return Response("No Centrifuge tasks. This is not the API you are looking for....", status=404)

    # If there is no MetaGenomicsMeta object return an error
    try:
        queryset = Metadata.objects.filter(task__id=job_master.id).last()
    except Metadata.DoesNotExist as e:
        return Response(e, status=404)

    if not queryset:
        return Response("No Metadata. This is not the task you are looking for....", status=404)

    number_of_reads = job_master.flowcell.number_reads

    reads_class = job_master.read_count
    # Percentage of reads classified
    percentage = round(reads_class / number_of_reads * 100, 2)
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

    return_list = [{"key": "Reads Sequenced: ", "value": number_of_reads},
                   {"key": "Reads Classified: ", "value": reads_class},
                   {"key": "Classified: ", "value": str(percentage) + "%"},
                   {"key": "Runtime: ", "value": runtime}]
    return Response(return_list)


@api_view(["GET"])
def centrifuge_sankey(request):
    # TODO refactor logic into centrifuge.py
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

    # The most up to dat task_id
    task_id = JobMaster.objects.filter(flowcell__id=flowcell_id, job_type__name="Metagenomics").order_by("id").last().id

    # ## Get the links for the sankey Diagram ###

    queryset = SankeyLink.objects.filter(task__id=task_id, barcode_name=selected_barcode).values()
    # If the queryset is empty, return an empty object
    if not queryset:
        return Response({}, status=204)
    # Create the dataframe from the results of the database query
    source_target_df = pd.DataFrame(list(queryset))
    # If the database is empty return an empty list
    if source_target_df.empty:
        return Response({}, status=204)
    # get a subset df of all the species rows
    temp_species_df = source_target_df[source_target_df["target_tax_level"] == "species"]
    # get species limit (default 50) of the largest species
    temp_species_df = temp_species_df.nlargest(species_limit, ["value"])
    # get the values for those species
    source_target_df = source_target_df[source_target_df["tax_id"].isin(temp_species_df["tax_id"])]
    source_target_df.set_index(["source", "target"], inplace=True)
    # Group by rows where the source and target are the same
    st_gb = source_target_df.groupby(["source", "target"])
    # Replace the value columns on dataframe with the sum of all the values of identical source target rows
    source_target_df["value"] = st_gb["value"].sum()

    source_target_df.reset_index(inplace=True)
    # Drop all duplicate rows, only need one new entry
    source_target_df.drop_duplicates(["source", "target"], keep="first", inplace=True)
    # Drop any rows where the source and the target are the same and don't keep them
    source_target_df = pd.concat([source_target_df, source_target_df[
        source_target_df["source"] == source_target_df["target"]]]).drop_duplicates(["source", "target"], keep=False)
    source_target_df.dropna(inplace=True)

    # Drop unnecessary columns from DataFrame
    source_target_df.drop(columns=["tax_id", "target_tax_level",
                                   "task_id"], inplace=True)
    # Set MultiIndex to group source to target
    source_target_df = pd.concat([source_target_df, source_target_df[
        source_target_df["source"] == source_target_df["target"]]]).drop_duplicates(["source", "target"], keep=False)
    # source_target_df.sort_values(["value"], ascending=False, inplace=True)
    # Create the links list of dicts to return
    links = source_target_df.to_dict(orient="records")

    # ## Get the nodes ###
    # Get all the nodes values from superkingdom (ex. Bacteria) to species ("E. Coli")
    # Create a series of all possible nodes
    nodes = source_target_df["source"].append(source_target_df["target"])

    # Remove duplicates
    nodes = pd.DataFrame({"name": nodes.unique()})

    merge = pd.merge(nodes, source_target_df, how="outer", right_on="source", left_on="name")

    merge = merge[["name", "path"]]

    merge = pd.merge(merge, source_target_df, how="outer", right_on="target", left_on="name")

    merge["path_x"].fillna(merge["path_y"], inplace=True)

    merge = merge[["name", "path_x"]]

    merge.drop_duplicates(subset="name", inplace=True, keep="first")

    merge.rename(columns={"path_x": "path"}, inplace=True)

    # Put the data in the right format
    nodes = merge.to_dict(orient="records")

    # ## Return the array # ##
    nodes = {"links": links, "nodes": nodes}
    return_dict = {"sankey": nodes}
    return Response(return_dict, status=200)


@api_view(["GET"])
def vis_table_or_donut_data(request):
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
    # Get the most recent job
    task = JobMaster.objects.filter(flowcell__id=flowcell_id, job_type__name="Metagenomics").order_by("id").last()
    results_df = pd.DataFrame()

    tax_rank_filter = ["superkingdom", "phylum", "classy", "order", "family", "genus", "species"]

    for rank in tax_rank_filter:
        data_df = pd.DataFrame(
            list(DonutData.objects.filter(task=task, tax_rank=rank,
                                          barcode_name=barcode).order_by("-num_matches")[:10].values()))
        results_df = results_df.append(data_df)

    if results_df.empty:
        return Response("No results", status=204)

    results_df = results_df[["tax_rank", "name", "num_matches"]]

    results_df = results_df.rename(columns={"name": "label", "num_matches": "value"})

    results_df.fillna("Unclassified", inplace=True)

    gb = results_df.groupby("tax_rank")

    return_dict = {}

    for name, group in gb:
        return_dict[name] = group.to_dict(orient="records")

    return Response(return_dict, status=200)


@api_view(["GET"])
def get_target_mapping(request):
    """
    Get the target species
    :param request:
    :return:
    """
    flowcell_id = request.GET.get("flowcellId", 0)

    barcode = request.GET.get("barcode", "All reads")

    if flowcell_id == 0:
        return Response("Flowcell id has failed to be delivered", status=404)

    # Get the most recent jobmaster id
    task_id = JobMaster.objects.filter(flowcell__id=flowcell_id, job_type__name="Metagenomics").order_by("id").last().id

    if barcode == "All reads":
        queryset = MappingResult.objects.filter(task__id=task_id, barcode_name=barcode).values()
        results_df = pd.DataFrame(list(queryset))
        if results_df.empty:
            return Response("No data has yet been produced for the target mappings", status=204)

    else:
        # Get a line for all the targets
        map_queryset = MappingResult.objects.filter(task__id=task_id).values()

        if not map_queryset:
            return Response("No data has yet been produced for the target mappings", status=204)
        # Get the Results for this barcode
        queryset = MappingResultsBarcoded.objects.filter(mapping_result__task__id=task_id,
                                                         barcode_name=barcode
                                                         ).values()
        # Create a dataframe for this barcode
        just_barcode_df = pd.DataFrame(list(queryset))
        if just_barcode_df.empty:
            return Response("No data has been found for this barcode. Perhaps an old barcode href is being queried",
                            status=404)
        # Barcode data only has species with results, need to add on any targets that have no results in that barcode

        all_targets_df = pd.DataFrame(list(map_queryset))
        all_targets_df[["num_mapped", "red_reads", "mapped_proportion_of_classified",
                        "red_reads_proportion_of_classified", "red_sum_unique",
                        "red_sum_unique_proportion_of_classified"]] = 0
        # results_df = pd.merge(all_targets_df, just_barcode_df, how="outer", left_on="tax_id", right_on="tax_id")
        results_df = just_barcode_df.append(all_targets_df)
        results_df.drop_duplicates(subset=["species"], keep="first", inplace=True)

    species_list = MappingResult.objects.filter(task__id=task_id).values_list("tax_id", flat=True)

    cent_output = CentrifugeOutput.objects.filter(task__id=task_id, tax_id__in=species_list,
                                                  barcode_name=barcode).values()
    if not cent_output:
        return Response("No data has been found for this barcode. Perhaps an old barcode href is being queried",
                        status=404)

    cent_output_df = pd.DataFrame(list(cent_output))

    if results_df.empty:
        queryset = MappingResult.objects.filter(task__id=task_id).values()

        results_df = pd.DataFrame(list(queryset))

        results_df[["num_mapped", "red_reads", "red_sum_unique"]] = 0

        results_df["Num. matches"] = 0

        results_df["Sum. Unique"] = 0

        results_df.rename(columns={"num_mapped": "Num. mapped",
                                   "red_reads": "Danger reads",
                                   "red_sum_unique": "Unique Danger reads",
                                   "species": "Species",
                                   "tax_id": "Tax id",
                                   "num_matches": "Num. matches",
                                   "sum_unique": "Sum. Unique",
                                   "mapped_proportion_of_classified": "Mapped prop. total",
                                   "red_reads_proportion_of_classified": "Red prop. total",
                                   "proportion_of_classified": "Prop. classified"
                                   }, inplace=True)

        results = results_df.to_dict(orient="records")

        return Response(results)

    merger_df = pd.merge(results_df, cent_output_df, how="outer", left_on="tax_id", right_on="tax_id")
    merger_df.drop(columns=["id_x", "id_y", "barcode_name_x", "barcode_name_y"], inplace=True)
    merger_df.fillna(0, inplace=True)

    merger_df.rename(columns={"num_mapped": "Num. mapped",
                              "red_reads": "Danger reads",
                              "red_sum_unique": "Unique Danger reads",
                              "species_x": "Species",
                              "tax_id": "Tax id",
                              "num_matches": "Num. matches",
                              "sum_unique": "Sum. Unique",
                              "mapped_proportion_of_classified": "Mapped prop. total (%)",
                              "red_reads_proportion_of_classified": "Red prop. total (%)",
                              "proportion_of_classified": "Prop. classified (%)"
                              }, inplace=True)

    results = merger_df.to_dict(orient="records")

    return Response(results)


@api_view(['GET'])
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

        task = JobMaster.objects.filter(flowcell=flowcell, job_type__name="Metagenomics").order_by('id').last()

        if task:
            metagenomics_barcodes = CentrifugeOutput.objects.filter(task__id=task.id) \
                .values_list("barcode_name", flat=True).distinct()

    return Response({"data": metagenomics_barcodes})


def all_results_table(request):
    """
    Returns the data for the metagenomics tab all results table
    :param request:
    :return:
    """
    query_columns = [
        'barcode_name',
        'superkingdom',
        'phylum',
        'classy',
        'order',
        'family',
        'genus',
        'species',
        'num_matches',
        'proportion_of_classified',
    ]

    flowcell_id = int(request.GET.get('flowcell_id', 5))

    draw = int(request.GET.get('draw', 0))

    start = int(request.GET.get('start', 0))

    length = int(request.GET.get('length', 10))

    end = start + length
    # Which column s
    search_value = request.GET.get('search[value]', '')

    order_column = request.GET.get('order[0][column]', '')
    # ascending descending
    order_dir = request.GET.get('order[0][dir]', '')

    if not search_value == "":

        cent_out_temp = CentrifugeOutput.objects \
            .filter(task__flowcell_id=flowcell_id) \
            .filter(task__flowcell__owner=request.user) \
            .filter(species__icontains=search_value) | CentrifugeOutput.objects \
            .filter(task__flowcell_id=flowcell_id) \
            .filter(task__flowcell__owner=request.user) \
            .filter(genus__icontains=search_value) | CentrifugeOutput.objects \
            .filter(task__flowcell_id=flowcell_id) \
            .filter(task__flowcell__owner=request.user) \
            .filter(family__icontains=search_value) | CentrifugeOutput.objects \
            .filter(task__flowcell_id=flowcell_id) \
            .filter(task__flowcell__owner=request.user) \
            .filter(order__icontains=search_value)

    else:

        cent_out_temp = CentrifugeOutput.objects \
            .filter(task__flowcell_id=flowcell_id) \
            .filter(task__flowcell__owner=request.user)

    if order_column:

        if order_dir == 'desc':
            cent_out_temp = cent_out_temp.order_by('-{}'.format(query_columns[int(order_column)]))

        else:
            cent_out_temp = cent_out_temp.order_by('{}'.format(query_columns[int(order_column)]))

    cents = cent_out_temp.values('tax_id',
                                 'barcode_name',
                                 'name',
                                 'num_matches',
                                 'proportion_of_classified',
                                 'superkingdom',
                                 'phylum',
                                 'classy',
                                 'order',
                                 'family',
                                 'genus',
                                 'species',
                                 )

    records_total = cent_out_temp.count()

    result = {
        'draw': draw,
        "recordsTotal": records_total,
        "recordsFiltered": records_total,
        "data": list(cents[start:end])
    }

    return JsonResponse(result, safe=True)
