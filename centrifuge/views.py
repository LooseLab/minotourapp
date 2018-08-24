from rest_framework.response import Response
from rest_framework.decorators import api_view
from centrifuge.models import CentOutput, CartographyMapped, CartographyGuide, LineageValues, MetaGenomicsMeta, \
    SankeyLinks
from centrifuge.serializers import CentSerialiser, CartMappedSerialiser, CartGuideSerialiser
from django.http import HttpResponse
from django.utils import timezone
import celery.bin.celery
import celery.platforms
import celery.bin.base
import uuid
from ete3 import NCBITaxa
import pandas as pd
from devices.models import Flowcell
from jobs.models import JobMaster, JobType
from datetime import datetime


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

    flowcell_id = request.GET.get("flowcellId", False)

    task_ids = JobMaster.objects.filter(flowcell=flowcell_id).values_list("id", flat=True)
    # Get the most recent job
    task_id = max(task_ids)
    try:
        queryset = MetaGenomicsMeta.objects.get(flowcell_id=flowcell_id,
                                            task__id=task_id)
    except MetaGenomicsMeta.DoesNotExist:
        return Response([], status=404)

    number_of_reads = queryset.number_of_reads
    reads_class = queryset.reads_classified
    # Percentage of reads classified
    percentage = reads_class / number_of_reads * 100
    start_time = queryset.timestamp.replace(tzinfo=None)
    # use the current time to work out how long the ongoing run has taken
    current_time = timezone.now().replace(tzinfo=None)
    # If the request is for a finished run, return the time taken using the set finish time
    if not queryset.finish_time:
        runtime = current_time - start_time
        runtime = str(runtime)
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
    Retrieves the data for the sankey diagram, returns it
    Parameters
    ----------
    request - The request body and parameters

    Returns - Response containing all the information for the the sankey diagram
    -------

    """
    species_limit = request.GET.get("speciesLimit", 50)
    flowcell_id = request.GET.get("flowcellId", "")
    task_id = max(JobMaster.objects.filter(flowcell=flowcell_id).values_list("id", flat=True))

    # ## Get the links for the sankey Diagram ###
    print(f"the flowcell is  {request.GET.get('flowcellId', '')}")
    queryset = SankeyLinks.objects.filter(flowcell_id=flowcell_id, task__id=task_id).values()

    if not queryset:
        return Response({}, status=204)

    source_target_df = pd.DataFrame(list(queryset))

    if source_target_df.empty:
        return Response({}, status=204)
    # Reduce the DataFrame to the number of species that you need
    source_target_df = source_target_df[source_target_df["rank"] < species_limit]
    # Drop unnecessary columns from DataFrame
    source_target_df.drop(columns=["flowcell_id", "id", "tax_id", "rank"], inplace=True)
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
        source_target_df["source"] == source_target_df["target"]]]).drop_duplicates(["source", "target"],
                                                                                    keep=False)
    source_target_df.dropna(inplace=True)
    links = source_target_df.to_dict(orient="records")

    # ## Get the nodes ###
    # Get all the nodes values from superkingdom (ex. Bacteria) to species ("E. Coli")
    nodes = source_target_df["source"].append(source_target_df["target"])

    # Remove duplicates
    nodes = pd.DataFrame({"name": nodes.unique()})

    # Put the data in the right format
    nodes = nodes.to_dict(orient="records")

    # ## Return the array # ##
    nodes = {"links": links, "nodes": nodes}
    return_dict = {"sankey": nodes}
    return Response(return_dict, status=200)


# TODO refactor into run centrifuge.py
@api_view(["GET"])
def vis_table_or_donut_data(request):
    """
    Create the total reads table for the bottom of the page
    Parameters
    ----------
    request - Django rest framework request object

    Returns - a dict containing a json string that has all the data for the total reads table
    -------

    """
    flowcell_id = request.GET.get("flowcellId", 0)
    # Get the most recent job
    task_id = max(JobMaster.objects.filter(flowcell=flowcell_id).values_list("id", flat=True))
    print(f"the flowcell job id is {task_id}, the flowcell id is {flowcell_id}")
    # queryset from database, filtered by the meta_id
    queryset = CentOutput.objects.filter(flowcell_id=flowcell_id, task__id=task_id)
    # all the taxIDs for this centrifuge job as a list, need for filtering all the relevant lineages
    centouput_df = pd.DataFrame(list(queryset.values()))
    # if there is no data in the database (yet) return 404
    if centouput_df.empty:
        return Response([], status=204)
    centouput_df.drop(columns=["task_id"], inplace=True)

    # create dataframe
    lineages_df = pd.DataFrame(list(LineageValues.objects.all().values()))
    # get relevant taxids, so dataframe of only taxIDs that are in the output of this centrifuge run
    lineages_df = lineages_df[lineages_df['tax_id'].isin(centouput_df["tax_id"])]
    # create the numreads series by mapping over the tax_id series
    lineages_df.set_index("tax_id", inplace=True)
    centouput_df.set_index("tax_id", inplace=True)
    lineages_df["num_matches"] = centouput_df["num_matches"]
    lineages_df["sum_unique"] = centouput_df["sum_unique"]

    # try and uniform missing values
    lineages_df = lineages_df.replace("nan", "Empty")
    lineages_df = lineages_df.replace("None", "Empty")
    # The order that the phyla go in
    # TODO add in substrain and subspecies when index contains them
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
        print(clade)
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

    lineages_df.fillna("Unknown", inplace=True)
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


@api_view(["POST"])
def start_centrifuge_view(request):
    """

    Parameters
    ----------
    request - (type DRF Request object) contains all the metadata we need to store about the run

    Returns - Response status 200 if celery task starts
    -------

    """
    # # create a new uuid, to use as a unique identifier
    # ident = str(uuid.uuid4())
    # # setup check to see if celery is up and running
    # status = celery.bin.celery.CeleryCommand.commands['status']()
    # status.app = status.get_app()
    # try:
    #     # see if celery is actually running
    #     status.run()
    #     # Create metadata
    #     data = request.data
    #     # get date and time
    #     d = datetime.now()
    #     # date is great
    #     # time is lime
    #     time = "{:%H:%M:%S}".format(d)
    #
    #     flowcell_id = request.data["flowcellId"]
    #
    #     # call celery task, pass the uuid we generated as an argument
    #     # Create Jobmaster object, to save. This will be picked up by the run
    #     # monitor task
    #     flowcell = Flowcell.objects.get(pk=flowcell_id)
    #     jobtype = JobType.objects.get(name="Centrifuge")
    #     task = JobMaster(flowcell=flowcell, job_type=jobtype,
    #                      running=False, complete=0, read_count=0, last_read=0)
    #     task.save()
    #     # create django model object, save it to database
    #     MetaGenomicsMeta(run_time=time, flowcell_id=flowcell_id, running=True, number_of_reads=0,
    #                      reads_classified=0, task=task).save()
    #     # # create django model object, save it to database
    #     # print("saved metadata")
    #     # return_dict = {"ident": ident}
    #     return Response(status=200)
    # # except if celery isn't running
    # except celery.bin.base.Error as e:
    #     if e.status == celery.platforms.EX_UNAVAILABLE:
    #         print("START CELERY")
    #         return HttpResponse(e).status_code(500)


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
