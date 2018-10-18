from rest_framework.response import Response
from rest_framework.decorators import api_view
from centrifuge.models import CentOutput, CartographyMapped, LineageValues, MetaGenomicsMeta, \
    SankeyLinks, CentOutputBarcoded, SankeyLinksBarcode
from centrifuge.serializers import CartMappedSerialiser
from django.utils import timezone
from ete3 import NCBITaxa
import pandas as pd
from jobs.models import JobMaster


@api_view(["GET"])
def metaview(request):
    """
    :purpose: return the metadata displayed at the top of the page relating to the task data eing currently visualised

    :author: Rory

    :param request: A django rest framework request body object, contains the flowcell id as a query parameter
    :return: A Response object, containg a list with the four data objects, Reads sequenced, Reads Classified,
    Classified and the Runtime
    """
    # The flowcell id for the flowcell that the fastq data came from, default to False if not present
    flowcell_id = request.GET.get("flowcellId", False)
    # The ids of the JobMaster entries for this flowcell
    task_ids = JobMaster.objects.filter(flowcell__id=flowcell_id, job_type__name="Metagenomics")\
        .values_list("id", flat=True)
    # Get the most recent job, which has the highest ID
    task_id = max(task_ids)
    # If there is no MetaGenomicsMeta object return an empty list
    try:
        queryset = MetaGenomicsMeta.objects.filter(flowcell__id=flowcell_id,
                                                task__id=task_id).last()
    except MetaGenomicsMeta.DoesNotExist:
        return Response([], status=404)

    number_of_reads = queryset.number_of_reads
    reads_class = queryset.reads_classified
    # Percentage of reads classified
    percentage = reads_class / number_of_reads * 100
    # Get the start time of the task, removing the timezone info
    start_time = queryset.timestamp.replace(tzinfo=None)
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
def cent_sankey_two(request):
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
    # Selected barcode, default all reads
    selected_barcode = request.GET.get("barcode", "All reads")

    # The most up to dat task_id
    task_id = max(JobMaster.objects.filter(flowcell__id=flowcell_id,
                                           job_type__name="Metagenomics").values_list("id", flat=True))

    # ## Get the links for the sankey Diagram ###

    print("the flowcell is {}".format(flowcell_id))
    print("the barcode is {}".format(selected_barcode))
    queryset = SankeyLinks.objects.filter(flowcell__id=flowcell_id, task__id=task_id, barcode=selected_barcode).values()
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

    # barcode_df = pd.DataFrame(list(SankeyLinksBarcode.objects.filter(link__flowcell_id=flowcell_id,
    #                                                                  link__task__id=task_id,
    #                                                                  barcode=selected_barcode).values()))
    # Merge barcode values onto dataframe
    # source_target_df = pd.merge(source_target_df, barcode_df, how="inner", on="tax_id")
    # get a subset df of all the species rows
    temp_species_df = source_target_df[source_target_df["target_tax_level"] == "species"]
    # get species limit (default 50) of the largest species
    temp_species_df = temp_species_df.nlargest(species_limit, ["value"])
    # get the values for those species
    source_target_df = source_target_df[source_target_df["tax_id"].isin(temp_species_df["tax_id"])]

    # Drop unnecessary columns from DataFrame
    source_target_df.drop(columns=["flowcell_id", "tax_id", "target_tax_level",
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


# TODO refactor into run centrifuge.py then save results in the correct format
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
    print()
    # Get the most recent job
    task_id = max(JobMaster.objects.filter(flowcell__id=flowcell_id, job_type__name="Metagenomics")
                  .values_list("id", flat=True))
    # queryset from database, filtered by the flowcell_id and the corresponding JobMaster ID

    queryset = CentOutputBarcoded.objects.filter(output__flowcell__id=flowcell_id, output__task__id=task_id,
                                                 barcode=barcode)
    # Create a dataframe from the results of the querying the database
    centoutput_df = pd.DataFrame(list(queryset.values()))
    # if there is no data in the database (yet) return 204
    if centoutput_df.empty:
        return Response([], status=204)
    names = CentOutput.objects.filter(flowcell__id=flowcell_id, task__id=task_id).values()

    names_df = pd.DataFrame(list(names))
    # if there is no names data
    if names_df.empty:
        return Response([], status=204)
    names_df.set_index("tax_id", inplace=True)
    centoutput_df.set_index("tax_id", inplace=True)

    centoutput_df["name"] = names_df["name"]
    centoutput_df.drop(columns=["output_id"], inplace=True)

    # create dataframe of all the Lineages we have stored in the database
    lineages_df = pd.DataFrame(list(LineageValues.objects.all().values()))
    # get relevant taxids, so dataframe of only lineages for taxIDs that are in the output of this centrifuge run
    lineages_df = lineages_df[lineages_df['tax_id'].isin(centoutput_df.index.values)]
    # Create the new columns in the dataframe for num_matches and num_unique matches
    lineages_df.set_index("tax_id", inplace=True)
    lineages_df["num_matches"] = centoutput_df["num_matches"]
    lineages_df["sum_unique"] = centoutput_df["sum_unique"]

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
        sumtitle = "summed_" + clade
        proptitle = "prop_" + clade
        # create new series, contains summed reads for this clade
        df[sumtitle] = gb
        # create new seies, contains proportion of reads in that clade across the whole dataframe
        df[proptitle] = df[sumtitle].div(df[sumtitle].unique().sum()).mul(100)
        # round to two decimal places
        df[proptitle] = df[proptitle].round(decimals=3)
        df = df.reset_index()
        return df

    def sort_top_ten_clade(clade, df):
        """
        sorts out the top 20 records in each clade, and transforms them into the right form for the donut chart
        Parameters
        ----------
        clade - (type String)
        df - (type pandas DataFrame)

        Returns - a list of dictionaries containg the top 20 members of a taxa format
        [{"label": "bacteria": "value": 36000}] where value is the num_reads in a clade
        -------

        """
        sumtitle = "summed_" + clade
        # create temp df with just the number of reads and the member of that clade
        temp_df = df[[clade, sumtitle]]
        temp_df.set_index(clade, inplace=True)
        temp_df.drop_duplicates(keep="first", inplace=True)
        temp_df.reset_index(inplace=True)
        # get the largest 20 members
        temp_df = temp_df.nlargest(10, sumtitle)
        # rename columns to what's expected by the donut chart
        temp_df.rename(columns={clade: "label", sumtitle: "value"}, inplace=True)
        # get the resulst in a list with adict for each record in temp_df
        list_dict_of_members = temp_df.to_dict(orient="records")
        return list_dict_of_members

    container_array = []

    # for each clade call the above function to perform the calculation and add the series
    for ord in order:
        lineages_df = sum_prop_clade_reads(ord, lineages_df)
        arr = sort_top_ten_clade(ord, lineages_df)
        obj = {ord: arr}
        container_array.append(obj)

    lineages_df.fillna("Unclassified", inplace=True)
    json = lineages_df.to_dict(orient="records")

    container_array.reverse()

    # create an object to return
    if request.GET.get("visType", "") == "table":
        return_dict = json
    elif request.GET.get("visType", "") == "donut":
        return_dict = {"result": container_array}
    else:
        return Response(status=400)

    # return it
    return Response(return_dict, status=200)


@api_view(["POST", "GET"])
def get_or_set_cartmap(request):
    """
    TODO Currently Unused
    Get the default mapping targets to display or set them for the run
    Parameters
    ----------
    request - (type DRF request object)

    Returns - either the list to be viewed in the start task table or a status that they have been inserted into db
    -------

    """
    if request.method == "GET":
        queryset = CartographyMapped.objects.filter(task_meta=request.query_params.get("meta_id", False))
        data = CartMappedSerialiser(queryset, many=True)
        return_data = []
        for data in data.data:
            data["amb_proportion"] = 0
            data["red_proportion"] = 0
            return_data.append(data)
        return_object = {"data": return_data}
        return Response(return_object, status=200)

    elif request.method == "POST":
        ncbi = NCBITaxa()

        for r in request.data["data"]:
            if type(r["species"]) is dict:
                del r
                continue
            ident = request.data["ident"]
            if "tax_id" not in r.keys():
                species = [r["species"]]
                name_to_taxid = ncbi.get_name_translator(species)
                tax_id = name_to_taxid[species[0]][0]
                CartographyMapped(species=species[0], tax_id=tax_id, task_meta=ident, total_reads=0, num_reads=0,
                                  alert_level=0, red_reads=0).save()
            else:
                CartographyMapped(species=r["species"], tax_id=r["tax_id"], task_meta=ident, total_reads=0, num_reads=0,
                                  alert_level=0, red_reads=0).save()

        return Response("Hi there", status=200)


@api_view(["POST", "GET"])
def get_target_mapping(request):
    flowcell_id = request.GET.get("flowcellId", 0)

    if flowcell_id == 0:
        return Response(status=404)

    print("flowcell_id_id-{}".format(flowcell_id))
    task_id = max(JobMaster.objects.filter(flowcell__id=flowcell_id, job_type__name="Metagenomics")
                  .values_list("id", flat=True))

    queryset = CartographyMapped.objects.filter(flowcell__id=flowcell_id,
                                                task__id=task_id).values()

    return Response(queryset)



