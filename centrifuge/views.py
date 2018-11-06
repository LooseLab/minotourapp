"""
views.py
"""
from rest_framework.response import Response
from rest_framework.decorators import api_view
from centrifuge.models import CentrifugeOutput, MappingResult, LineageValue, Metadata, \
    SankeyLink, CentrifugeOutputBarcoded, MappingResultsBarcoded
from django.utils import timezone

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
        return Response("No Centrifuge tasks. This is not the API you are looking for....")

    # If there is no MetaGenomicsMeta object return an error
    try:
        queryset = Metadata.objects.filter(task__id=job_master.id).last()
    except Metadata.DoesNotExist as e:
        return Response(e, status=404)
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
    species_limit = request.GET.get("speciesLimit", 30)
    # Get the flowcell ID , defaulting to False if not present
    flowcell_id = request.GET.get("flowcellId", False)
    print(flowcell_id)
    # Selected barcode, default all reads
    selected_barcode = request.GET.get("barcode", "All reads")

    # The most up to dat task_id
    task_id = JobMaster.objects.filter(flowcell__id=flowcell_id, job_type__name="Metagenomics").order_by("id").last().id

    # ## Get the links for the sankey Diagram ###

    print("the flowcell is {}".format(flowcell_id))
    print("the barcode is {}".format(selected_barcode))
    queryset = SankeyLink.objects.filter(task__id=task_id, barcode_name=selected_barcode).values()
    # If the queryset is empty, return an empty object
    if not queryset:
        print("no queryset")
        return Response({}, status=204)
    # Create the dataframe from the results of the database query
    source_target_df = pd.DataFrame(list(queryset))
    # If the database is empty return an empty list
    if source_target_df.empty:
        print("no source target")
        return Response({}, status=204)
    # get a subset df of all the species rows
    temp_species_df = source_target_df[source_target_df["target_tax_level"] == "species"]
    # get species limit (default 50) of the largest species
    temp_species_df = temp_species_df.nlargest(species_limit, ["value"])
    # get the values for those species
    source_target_df = source_target_df[source_target_df["tax_id"].isin(temp_species_df["tax_id"])]

    # Drop unnecessary columns from DataFrame
    source_target_df.drop(columns=["tax_id", "target_tax_level",
                                   "task_id"], inplace=True)
    # Set MultiIndex to group source to target
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
    source_target_df.sort_values(["value"], ascending=False, inplace=True)
    # Create the links list of dicts to return
    links = source_target_df.to_dict(orient="records")

    # ## Get the nodes ###
    # Get all the nodes values from superkingdom (ex. Bacteria) to species ("E. Coli")
    # Create a series of all possible nodes
    nodes = source_target_df["source"].append(source_target_df["target"])

    # Remove duplicates
    nodes = pd.DataFrame({"name": nodes.unique()})

    # Put the data in the right format
    nodes = nodes.to_dict(orient="records")

    # ## Return the array # ##
    nodes = {"links": links, "nodes": nodes}
    return_dict = {"sankey": nodes}
    return Response(return_dict, status=200)


# TODO refactor into run centrifuge.py then save results in the correct format, cuase this is a slow point
# TODO refactor table into server side preprocessing
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
    task_id = JobMaster.objects.filter(flowcell__id=flowcell_id, job_type__name="Metagenomics").order_by("id").last().id
    # queryset from database, filtered by the flowcell_id and the corresponding JobMaster ID

    queryset = CentrifugeOutputBarcoded.objects.filter(
        output__task__id=task_id,
        barcode_name=barcode)

    barcode_list = list(CentrifugeOutputBarcoded.objects.filter(
        output__task__id=task_id)
                        .values_list("barcode_name", flat=True).distinct())

    # Create a dataframe from the results of the querying the database
    centrifuge_output_barcoded_df = pd.DataFrame(list(queryset.values()))

    # if there is no data in the database (yet) return 204
    if centrifuge_output_barcoded_df.empty:
        return Response([], status=204)

    species_names = CentrifugeOutput.objects.filter(task__id=task_id).values()

    species_names_df = pd.DataFrame(list(species_names))
    # if there is no species_names data
    if species_names_df.empty:
        return Response([], status=204)

    species_names_df.set_index("tax_id", inplace=True)

    centrifuge_output_barcoded_df.set_index("tax_id", inplace=True)

    centrifuge_output_barcoded_df["name"] = species_names_df["name"]

    centrifuge_output_barcoded_df.drop(columns=["output_id"], inplace=True)

    # create dataframe of all the Lineages we have stored in the database
    lineages_df = pd.DataFrame(list(LineageValue.objects.all().values()))
    # get relevant taxids, so dataframe of only lineages for taxIDs that are in the output of this centrifuge run
    lineages_df = lineages_df[lineages_df['tax_id'].isin(centrifuge_output_barcoded_df.index.values)]
    # Create the new columns in the dataframe for num_matches and num_unique matches
    lineages_df.set_index("tax_id", inplace=True)

    lineages_df["num_matches"] = centrifuge_output_barcoded_df["num_matches"]

    lineages_df["sum_unique"] = centrifuge_output_barcoded_df["sum_unique"]

    # try and uniform missing values
    lineages_df = lineages_df.replace("nan", "Unclassified")

    lineages_df = lineages_df.replace("None", "Unclassified")

    lineages_df = lineages_df.fillna("Unclassified")
    # The order that the phyla go in
    # TODO add in substrain and subspecies when centrifuge index contains them
    order = ["superkingdom", "phylum", "classy", "order", "family", "genus", "species"]

    def sum_prop_clade_reads(clade, df):
        """
        I wrote this function in Manchester! It rained. It calculates the sum and proportion of reads within a clade,
        and returns them into a new series
        Parameters
        ----------
        clade - (type String) - A taxonomic clade, i.e kingdom, phylum etc.
        df - (type pandas DataFrame) - The dataframe containing all the centrifuge output objects

        Returns
        -------
        The new Dataframe
        """
        # set the index so we can map the group by object back
        df = df.set_index(clade)
        # group the records by clade
        gb = df.groupby(level=clade)
        # sum number of reads in each clade
        gb = gb["num_matches"].sum()

        sum_title = "summed_" + clade

        proportion_title = "prop_" + clade
        # create new series, contains summed reads for this clade
        df[sum_title] = gb
        # create new seies, contains proportion of reads in that clade across the whole dataframe
        df[proportion_title] = df[sum_title].div(df[sum_title].unique().sum()).mul(100)
        # round to two decimal places
        df[proportion_title] = df[proportion_title].round(decimals=3)
        df = df.reset_index()
        return df

    def sort_top_ten_clade(clade, df):
        """
        sorts out the top 10 records in each clade, and transforms them into the right form for the donut chart
        Parameters
        ----------
        clade - (type String)
        df - (type pandas DataFrame)

        Returns - a list of dictionaries containg the top 20 members of a taxa format
        [{"label": "bacteria": "value": 36000}] where value is the num_reads in a clade
        -------

        """
        summed_title = "summed_" + clade
        # create temp df with just the number of reads and the member of that clade
        temp_df = df[[clade, summed_title]]

        temp_df.set_index(clade, inplace=True)

        temp_df.drop_duplicates(keep="first", inplace=True)

        temp_df.reset_index(inplace=True)
        # get the largest 20 members
        temp_df = temp_df.nlargest(10, summed_title)
        # rename columns to what's expected by the donut chart
        temp_df.rename(columns={clade: "label", summed_title: "value"}, inplace=True)
        # get the resulst in a list with adict for each record in temp_df
        list_dict_of_top_ten_members = temp_df.to_dict(orient="records")

        return list_dict_of_top_ten_members

    container_array = []

    # for each clade call the above function to perform the calculation and add the series
    for ordy in order:
        lineages_df = sum_prop_clade_reads(ordy, lineages_df)

        arr = sort_top_ten_clade(ordy, lineages_df)

        obj = {ordy: arr}

        container_array.append(obj)

    lineages_df.fillna("Unclassified", inplace=True)

    json = lineages_df.to_dict(orient="records")
    # Species to kingdom
    container_array.reverse()

    # create an object to return
    if request.GET.get("visType", "") == "table":
        return_dict = json

    elif request.GET.get("visType", "") == "donut":
        return_dict = {"result": container_array, "barcodes": barcode_list}

    else:
        return Response(status=400)

    # return it
    return Response(return_dict, status=200)


@api_view(["POST", "GET"])
def get_target_mapping(request):
    """
    Get the target species
    :param request:
    :return:
    """
    flowcell_id = request.GET.get("flowcellId", 0)

    barcode = request.GET.get("barcode", "All reads")

    if flowcell_id == 0:
        return Response(status=404)

    print("flowcell_id_id-{}".format(flowcell_id))

    task_id = JobMaster.objects.filter(flowcell__id=flowcell_id, job_type__name="Metagenomics").order_by("id").last().id

    if barcode == "All reads":
        queryset = MappingResult.objects.filter(task__id=task_id, barcode_name=barcode).values()

        results_df = pd.DataFrame(list(queryset))
    else:
        map_queryset = MappingResult.objects.filter(task__id=task_id).values()

        queryset = MappingResultsBarcoded.objects.filter(mapping_result__task__id=task_id,
                                                         barcode_name=barcode
                                                         ).values()
        just_barcode_df = pd.DataFrame(list(queryset))
        # Barcode data only has species with results, need to add on any targets that have no results in that barcode

        all_targets_df = pd.DataFrame(list(map_queryset))
        print(all_targets_df)
        print(just_barcode_df)
        # results_df = pd.merge(all_targets_df, just_barcode_df, how="outer", left_on="tax_id", right_on="tax_id")
        results_df = just_barcode_df.append(all_targets_df)
        results_df.drop_duplicates(subset=["species"], keep="first", inplace=True)
        print(results_df)

    species_list = MappingResult.objects.filter(task__id=task_id).values_list("tax_id", flat=True)

    cent_output = CentrifugeOutputBarcoded.objects.filter(output__task__id=task_id, tax_id__in=species_list,
                                                          barcode_name="All reads").values()

    print(results_df)
    cent_output_df = pd.DataFrame(list(cent_output))
    # TODO sorted
    print(cent_output_df)
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
                                   "sum_unique": "Sum. Unique"
                                   }, inplace=True)

        results = results_df.to_dict(orient="records")
        print(results_df)

        return Response(results)

    merger_df = pd.merge(results_df, cent_output_df, how="outer", left_on="tax_id", right_on="tax_id")
    merger_df.drop(columns=["id_x", "id_y", "barcode_name_x", "barcode_name_y", "output_id"], inplace=True)
    merger_df.fillna(0, inplace=True)
    print(merger_df)

    merger_df.rename(columns={"num_mapped": "Num. mapped",
                              "red_reads": "Danger reads",
                              "red_sum_unique": "Unique Danger reads",
                              "species": "Species",
                              "tax_id": "Tax id",
                              "num_matches": "Num. matches",
                              "sum_unique": "Sum. Unique"
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
    print("Metagenomics barcodes")
    flowcell_list = Flowcell.objects.filter(owner=request.user).filter(id=pk)
    print(pk)
    metagenomics_barcodes = []

    if flowcell_list.count() > 0:

        flowcell = flowcell_list[0]

        task = JobMaster.objects.filter(flowcell=flowcell, job_type__name="Metagenomics").order_by('id').last()

        if task:
            metagenomics_barcodes = CentrifugeOutputBarcoded.objects.filter(output__task__id=task.id)\
                .values_list("barcode_name", flat=True).distinct()

    return Response({"data": metagenomics_barcodes})
