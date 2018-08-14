from rest_framework import viewsets
from rest_framework.response import Response
from rest_framework.decorators import api_view
from centrifuge.models import CentOutput, CartographyMapped, CartographyGuide, ReferenceGenomes, \
    LineageValues, MetaGenomicsMeta
from centrifuge.serializers import CentSerialiser, CartMappedSerialiser, CartGuideSerialiser, \
    ReferenceGenomeSerialiser
from django.views.decorators.csrf import csrf_exempt
from django.http import HttpResponse
from django.utils import timezone
import celery.bin.celery
import celery.platforms
import celery.bin.base
from operator import itemgetter
import uuid
from ete3 import NCBITaxa
import pandas as pd
from devices.models import Flowcell
from jobs.models import JobMaster, JobType
from datetime import datetime

class CentViewSet(viewsets.ModelViewSet):
    """Viewset for viewing and editing centrifuge output objects"""

    queryset = CentOutput.objects.all()
    serializer_class = CentSerialiser


class DefaultViewSet(viewsets.ModelViewSet):
    """Viewset for viewing and editing the default targets in the traffic light analysis"""

    queryset = CartographyGuide.objects.all()
    serializer_class = CartGuideSerialiser


class RefViewSet(viewsets.ModelViewSet):
    # TODO this will need rewriting so people can upload private refences and retrieve them by their userID

    queryset = ReferenceGenomes.objects.all()
    serializer_class = ReferenceGenomeSerialiser


@api_view(["GET"])
def metaview(request):
    """
    Returns the metadata relating to a run, for the header on the vis page
    Parameters
    ----------
    request - django rest framework request object from AJAX get

    Returns - A response containing the Run_time, amount of reads in the database, and how many have been classfied, and a percentage
    -------

    """
    print(f"meta_id is {request.query_params.get('flowcellId', False)}")
    # print(request.GET.flowcell_id)
    # the relevant metagenomics database entry
    queryset = MetaGenomicsMeta.objects.get(flowcell_id=request.query_params.get("flowcellId", False))
    # If the run has finished set the finish time permanently
    if not queryset.running and not queryset.finish_time:
        queryset.finish_time = timezone.now().replace(tzinfo=None)
        queryset.save()

    number_of_reads = queryset.number_of_reads
    reads_class = queryset.reads_classified
    # Percentage of reads classified
    percentage = reads_class / number_of_reads * 100
    start_time = queryset.timestamp.replace(tzinfo=None)
    # use the current time to work out how long the ongoing run has taken
    current_time = timezone.now().replace(tzinfo=None)
    # If the request is for a finished run, return the time taken using the set finish time
    if not queryset.running:
        runtime = queryset.finish_time.replace(tzinfo=None) - start_time
        runtime = str(runtime)
    # Else return the time taken using the current time
    else:
        runtime = current_time - start_time
        runtime = str(runtime)

    return_list = [{"key": "Reads Sequenced: ", "value": number_of_reads},
                   {"key": "Reads Classified: ", "value": reads_class},
                   {"key": "Classified: ", "value": str(percentage) + "%"},
                   {"key": "Runtime: ", "value": runtime}]
    return Response(return_list)


#TODO REFACTOR into one function - very repetitive
@api_view(["GET"])
def cent_sankey(request):
    """
    Calculates the data for the sankey diagram, returns it
    Parameters
    ----------
    request - The request body and parameters

    Returns - Response containing all the information for the table and the sankey diagram
    -------

    """

    # queryset from database, filtered by the meta_id
    queryset = CentOutput.objects.filter(task_meta=request.query_params.get("flowcellId", False))
    # all the taxIDs for this centrifuge job as a list, need for filtering all the relevant lineages
    cent_output_tax_ids = queryset.values_list(
        "taxID", flat=True)
    # if there is no data in the database (yet) return 404 # TODO better status code
    if len(queryset) == 0:
        return Response(status=404)
    # Deserialise Django model object
    serialiser = CentSerialiser(queryset, many=True)
    # the data from desrialisation
    serialised_data = serialiser.data
    # get all the lineages
    lineages = list(LineageValues.objects.all().values())
    # lookup dict so you can get the data corrsponding to a specified taxID
    lookup_dict = {x["taxID"]: x for x in serialised_data}

    def return_reads(tax_id):
        """

        Parameters
        ----------
        tax_id (type integer) - the taxID from the dataframe series

        Returns (type integer) - the number of reads from the lookup dict from that taxID
        -------

        """
        if tax_id in lookup_dict:
            return lookup_dict[tax_id]["num_reads"]

    def return_unique_reads(tax_id):
        """

        Parameters
        ----------
        tax_id (type integer) - the taxID from the dataframe series

        Returns (type integer) - the number of unique reads from the lookup dict for that taxID
        -------

        """
        if tax_id in lookup_dict:
            return lookup_dict[tax_id]["sum_unique"]

    # Create Nodes dictionary, contains nodes and links array, needed for d3 sankey
    # nodes key is a list of dicts withe format {name: node_name}
    # links key is a list of dicts with format {source: starting_node, target: finishing node, value: num_reads}
    nodes = {"nodes": [], "links": []}
    # Taxas, in the order that links must flow
    order = ["superkingdom", "phylum", "classy", "order", "family", "genus", "species", "subspecies", "strain"]
    # Create a dataframe from the lineages list result from the database
    df = pd.DataFrame(lineages)
    # get relevant taxids, so dataframe of only taxIDs that are in the output of this centrifuge run
    df = df[df['tax_id'].isin(cent_output_tax_ids)]
    # create the numreads series by mapping over the tax_id series
    df["num_reads"] = df["tax_id"].map(return_reads, na_action=None)
    # create the sumunqiue series by mapping over thetax_id series
    df["sum_unique"] = df["tax_id"].map(return_unique_reads)
    # try and uniform missing values
    df = df.replace("nan", "Empty")
    df = df.replace("None", "Empty")
    # iterate over the whole data frame
    # TODO rewrite as the other two to avoid unnecessary iteration
    for ind, row in df.iterrows():
        ordered = []
        # get the flow that skips any missing values
        for orde in order:
            if row[orde] != "Empty" and row[orde] != None:
                ordered.append(orde)
        # enumerate over the ordered list so links flow from one to the next
        for index, ori in enumerate(ordered):
            if row[ori] != "Empty" and row[ori] != None:
                # the beginning of the link
                source_lineage_value = row[ori]
                # the node name isn't already in the list, add it
                if {"name": source_lineage_value} not in nodes["nodes"]:
                    nodes["nodes"].append({"name": source_lineage_value})

                # if the index for the target won't be out of range
                if index != len(ordered) - 1:
                    target_taxa_index = index + 1
                    # if the target index exists
                    if row[ordered[target_taxa_index]] != "Empty" and row[ordered[target_taxa_index]] != None:
                        # the target node value
                        target_lineage_value = row[ordered[target_taxa_index]]
                        # link_id is the unique id used to find the index
                        link_id = source_lineage_value + target_lineage_value
                        # create the link dict
                        link = {"source": source_lineage_value, "target": target_lineage_value,
                                "value": row["num_reads"], "id": link_id,
                                "target_tax_level": ordered[target_taxa_index]}
                        # check to see if the link is already in the links list
                        index = next((index for (index, d) in enumerate(nodes["links"]) if d["id"]
                                      == link_id), None)
                        # if the index is not already in the list
                        if index == None:
                            # if index:
                            index = -1
                        # if link is already in the list, add the value of the new link to the existing link
                        if index > -1:
                            nodes["links"][index]["value"] += link["value"]
                        # if link isn't already in, add it in
                        else:
                            nodes["links"].append(link)
                            # create the values of the reads

    sotted = sorted(nodes["links"], key=itemgetter('value'), reverse=True)

    # create an object to return
    return_dict = {"queryset": nodes, "values": sotted, "serialiser_data": serialised_data}
    # return it
    return Response(return_dict, status=200)


@api_view(["GET"])
def vis_table_data(request):
    """
    Create the total reads table for the bottom of the page
    Parameters
    ----------
    request - Django rest framework request object

    Returns - a dict containing a json string that has all the data for the total reads table
    -------

    """
    # queryset from database, filtered by the meta_id
    queryset = CentOutput.objects.filter(task_meta=request.query_params.get("flowcellId", False))
    # all the taxIDs for this centrifuge job as a list, need for filtering all the relevant lineages
    cent_output_tax_ids = queryset.values_list(
        "taxID", flat=True)
    # if there is no data in the database (yet) return 404 # TODO better status code
    if len(queryset) == 0:
        return Response(status=404)
    # Deserialise Django model object
    serialiser = CentSerialiser(queryset, many=True)
    # the data from deserialisation
    serialised_data = serialiser.data
    # get all the lineages
    lineages = list(LineageValues.objects.all().values())
    # lookup dict so you can get the data corrsponding to a specified taxID
    lookup_dict = {x["taxID"]: x for x in serialised_data}

    def return_reads(tax_id):
        """

        Parameters
        ----------
        tax_id (type integer) - the taxID from the dataframe series

        Returns (type integer) - the number of reads from the lookup dict from that taxID
        -------

        """
        if tax_id in lookup_dict:
            return lookup_dict[tax_id]["num_reads"]

    def return_unique_reads(tax_id):
        """

        Parameters
        ----------
        tax_id (type integer) - the taxID from the dataframe series

        Returns (type integer) - the number of unique reads from the lookup dict for that taxID
        -------

        """
        if tax_id in lookup_dict:
            return lookup_dict[tax_id]["sum_unique"]

    # create dataframe
    df = pd.DataFrame(lineages)
    # get relevant taxids, so dataframe of only taxIDs that are in the output of this centrifuge run
    df = df[df['tax_id'].isin(cent_output_tax_ids)]
    # create the numreads series by mapping over the tax_id series
    df["num_reads"] = df["tax_id"].map(return_reads, na_action=None)
    # create the sumunqiue series by mapping over thetax_id series
    df["sum_unique"] = df["tax_id"].map(return_unique_reads)
    # try and uniform missing values
    df = df.replace("nan", "Empty")
    df = df.replace("None", "Empty")
    # The order that the phyla go in
    order = ["superkingdom", "phylum", "classy", "order", "family", "genus", "species", "subspecies", "strain"]

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

        """
        # setthe index so we can map the group by object back
        df = df.set_index(clade)
        # group the records by clade
        gb = df.groupby(level=clade)
        # sum number of reads in each clade
        gb = gb["num_reads"].sum()
        sumtitle = "summed_" + clade
        proptitle = "prop_" + clade
        # create new series, contains summed reads for this clade
        df[sumtitle] = gb
        # create new seies, contains proportion of reads in that clade across the whole dataframe
        df[proptitle] = df[sumtitle].div(df[sumtitle].unique().sum()).mul(100)
        # round to two decimal places
        df[proptitle] = df[proptitle].round(decimals=2)
        df = df.reset_index()
        return df
    # for each clade call the above function to perform the calculation and add the series
    for ord in order:
        df = sum_prop_clade_reads(ord, df)
    df["species"] = df["species"] + " (Num Reads - " + df["summed_species"].map(str) \
                    + " Proportion " + df["prop_species"].map(str) + "%)"
    df["genus"] = df["genus"] + " (Num Reads - " + df["summed_genus"].map(str) + " Proportion " \
                  + df["prop_genus"].map(str) + "%)"
    json = df.to_json(orient="records")
    # create an object to return
    return_dict = {"json": json}
    # return it
    return Response(return_dict, status=200)


@api_view(["GET"])
def cent_donut(request):
    """
    Create the donut chart data and return it
    Parameters
    ----------
    request - (type DRF request)

    Returns - dict containing an array of objects, one for each clade
    -------

    """
    # get the relevant centOutput database entries for this analysis
    queryset = CentOutput.objects.filter(task_meta=request.query_params.get("flowcellId", False))
    # if there is no data return a 401 and a message
    if not queryset:
        return Response({"message": "no data"}, status=401)
    # get all the tax_ids from the database query as a list
    cent_output_tax_ids = queryset.values_list(
        "taxID", flat=True)
    # serialise the database entries into a python object
    serialiser = CentSerialiser(queryset, many=True)
    # get the serialisation result data
    serialised_data = serialiser.data
    # get all the Lineages we have already calculated as a list
    # TODO ineffecient - memory expensive. Look at only using lineages we need, same in all 3,
    # TODO a lot of not DRY code here
    lineages = list(LineageValues.objects.all().values())
    # lookup dict so you can get the data corrsponding to a specified taxID
    lookup_dict = {x["taxID"]: x for x in serialised_data}

    def return_reads(tax_id):
        """

        Parameters
        ----------
        tax_id (type integer) - the taxID from the dataframe series

        Returns (type integer) - the number of reads from the lookup dict from that taxID
        -------

        """
        if tax_id in lookup_dict:
            return lookup_dict[tax_id]["num_reads"]

    def return_unique_reads(tax_id):
        """

        Parameters
        ----------
        tax_id (type integer) - the taxID from the dataframe series

        Returns (type integer) - the number of unique reads from the lookup dict for that taxID
        -------

        """
        if tax_id in lookup_dict:
            return lookup_dict[tax_id]["sum_unique"]

    # create dataframe
    df = pd.DataFrame(lineages)
    # get relevant taxids, so dataframe of only taxIDs that are in the output of this centrifuge run
    df = df[df['tax_id'].isin(cent_output_tax_ids)]
    # create the numreads series by mapping over the tax_id series
    df["num_reads"] = df["tax_id"].map(return_reads, na_action=None)
    # create the sumunqiue series by mapping over thetax_id series
    df["sum_unique"] = df["tax_id"].map(return_unique_reads)
    # try and uniform missing values # TODo does this even work
    df = df.replace("nan", "Empty")
    df = df.replace("None", "Empty")
    df = df.fillna("Empty")
    # nodes = {"nodes": [], "links": []}
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

        """
        # setthe index so we can map the group by object back
        df = df.set_index(clade)
        # group the records by clade
        gb = df.groupby(level=clade)
        # sum number of reads in each clade
        gb = gb["num_reads"].sum()
        sumtitle = "summed_" + clade
        proptitle = "prop_" + clade
        # create new series, contains summed reads for this clade
        df[sumtitle] = gb
        # create new seies, contains proportion of reads in that clade across the whole dataframe
        df[proptitle] = df[sumtitle].div(df[sumtitle].unique().sum()).mul(100)
        # round to two decimal places
        df[proptitle] = df[proptitle].round(decimals=2)
        df = df.reset_index()
        return df
    # container array for results to be sent back
    container_array = []

    def sort_top_twenty_clade(clade, df):
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
        print(clade)
        sumtitle = "summed_" + clade
        # create temp df with just the number of reeads and the member of that clade
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
    # for each clade run the two functions above
    for ord in order:
        df = sum_prop_clade_reads(ord, df)
        arr = sort_top_twenty_clade(ord, df)
        obj = {ord: arr}
        container_array.append(obj)
    # reverse the array so species appears first
    container_array.reverse()

    result_dict = {"result": container_array}
    # return_dict = {"counted": countable_result_cont}
    return Response(result_dict, status=200)


# stop csrf blocking by django, # TODO fix so it's done properly with authentication
@csrf_exempt
@api_view(["POST"])
def start_centrifuge_view(request):
    """

    Parameters
    ----------
    request - (type DRF Request object) contains all the metadata we need to store about the run

    Returns - Response status 200 if celery task starts
    -------

    """
    # create a new uuid, to use as a unique identifier
    ident = str(uuid.uuid4())
    # setup check to see if celery is up and running
    status = celery.bin.celery.CeleryCommand.commands['status']()
    status.app = status.get_app()
    try:
        # see if celery is actually running
        status.run()
        # Create metadata
        data = request.data
        # get date and time
        d = datetime.now()
        # date is great
        # time is lime
        time = "{:%H:%M:%S}".format(d)

        flowcell_id = request.data["flowcellId"]
        # create django model object, save it to database
        MetaGenomicsMeta(run_time=time, flowcell_id=flowcell_id ,running=True, number_of_reads=0,
                         reads_classified=0).save()
        # call celery task, pass the uuid we generated as an argument
        # Create Jobmaster object, to save. This will be picked up by the run
        # monitor task
        flowcell = Flowcell.objects.get(pk=flowcell_id)
        jobtype = JobType.objects.get(pk=10)
        JobMaster(flowcell=flowcell, job_type=jobtype,
                  running=False, complete=0, read_count=0, last_read=0).save()
        # # create django model object, save it to database
        # print("saved metadata")
        # return_dict = {"ident": ident}
        return Response(status=200)
    # except if celery isn't running
    except celery.bin.base.Error as e:
        if e.status == celery.platforms.EX_UNAVAILABLE:
            print("START CELERY")
            return HttpResponse(e).status_code(500)


@csrf_exempt
@api_view(["POST", "GET"])
def get_or_set_cartmap(request):
    """
    Currently Unused
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
