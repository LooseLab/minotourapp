import datetime
import json
from collections import defaultdict
from pprint import pprint

import numpy as np
import pandas as pd
from celery.utils.log import get_task_logger
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.core.paginator import EmptyPage, PageNotAnInteger, Paginator
from django.core.serializers.json import DjangoJSONEncoder
from django.db.models import Max, Q
from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404
from django.shortcuts import render
from django.utils import timezone
from guardian.shortcuts import assign_perm, remove_perm
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from alignment.models import PafRoughCov
from artic.models import ArticBarcodeMetadata
from centrifuge.models import CentrifugeOutput
from minotourapp import settings
from reads.models import (
    Barcode,
    FastqFile,
    FastqRead,
    FastqReadType,
    Flowcell,
    FlowcellChannelSummary,
    FlowcellHistogramSummary,
    FlowcellStatisticBarcode,
    FlowcellSummaryBarcode,
    GroupRun,
    Minion,
    MinionControl,
    MinionEvent,
    MinionEventType,
    MinionMessage,
    MinionRunStats,
    MinionRunInfo,
    MinionScripts,
    MinionInfo,
    Run,
    JobType,
    JobMaster,
)
from reads.serializers import (
    BarcodeSerializer,
    ChannelSummary,
    ChannelSummarySerializer,
    FastqFileSerializer,
    FastqReadGetSerializer,
    FastqReadTypeSerializer,
    FlowcellSerializer,
    FlowcellSummaryBarcodeSerializer,
    GroupRunSerializer,
    MinIONControlSerializer,
    MinIONEventSerializer,
    MinIONEventTypeSerializer,
    MinionMessageSerializer,
    MinionRunStatsSerializer,
    MinionRunInfoSerializer,
    MinIONScriptsSerializer,
    MinionSerializer,
    MinionInfoSerializer,
    RunSerializer,
    RunStatisticBarcode,
    RunStatisticBarcodeSerializer,
    RunSummaryBarcode,
    RunSummaryBarcodeSerializer,
    JobMasterSerializer,
    JobMasterInsertSerializer,
)
from reads.tasks.delete_flowcell_task import clear_artic_data
from reads.tasks.redis_tasks_functions import save_reads_bulk
from reads.utils import get_coords, return_temp_empty_summary, pause_job, clear_artic_command_job_masters
from readuntil.models import ExpectedBenefitChromosomes
from reference.models import ReferenceInfo
from web.delete_tasks import (
    delete_metagenomics_task,
    delete_alignment_task,
)

logger = get_task_logger(__name__)


@api_view(["GET"])
def proportion_of_total_reads_in_barcode_list(request, pk):
    """
    Return a list of the proportion of the total make up of a flowcells read data each barcode is.

    Author
    ------
    Adonis

    Parameters
    ----------
    request: rest_framework.request.Request
        Django rest framework object.
    pk: str
        Primary key of the flowcell
    Returns`
    -------
    list of dicts
        The name of the series, and the series value..

    """

    if not pk:
        return Response("No Flowcell primary key provided", status=status.HTTP_400_BAD_REQUEST)
    flowcell = Flowcell.objects.get(pk=pk)
    if flowcell.number_barcodes > 3:
        qs = FlowcellSummaryBarcode.objects.filter(flowcell=flowcell).exclude(
            barcode_name__in=["No barcode", "All reads"]
        )
        df = pd.DataFrame(
            qs.values(
                "barcode_name",
                "read_count",
                "rejection_status",
                "read_type_name",
                "status",
            )
        )
        df["status"] = np.where(df["status"] == "True", "Pass", "Fail")
        total_read_count = df["read_count"].sum()
        df["data"] = round(df["read_count"] / total_read_count * 100, 2)
        df["name"] = df["barcode_name"]
        categories = df["name"].unique().tolist()
        # get the barcodes that have a pass status and no corresponding fail data
        no_fail_pass_df = df[
            (df["status"] == "Pass") & (~df["barcode_name"].isin(df[df["status"] == "Fail"]["barcode_name"].values))]
        no_fail_pass_df["status"] = "Fail"
        no_fail_pass_df["data"] = 0
        # append new entries with 0
        # Get the barcode rows that have a fail status and no cprreponding pass data
        df = df.append(no_fail_pass_df)
        no_pass_fail_df = df[
            (df["status"] == "Fail") & (~df["barcode_name"].isin(df[df["status"] == "Pass"]["barcode_name"].values))]
        no_pass_fail_df["status"] = "Pass"
        no_pass_fail_df["data"] = 0
        df = df.append(no_pass_fail_df)
        df = df.sort_values("barcode_name")
        # Get the proportion unclassified for pie chart
        df["is_unclassified"] = np.where(df["barcode_name"] == "unclassified", "Classified", "Unclassified")
        # use this array for easy comparison of current data on client side in BaseCalledData controller
        array_data = []
        pie_chart_data = []
        for key, value in df.groupby("is_unclassified")["data"].sum().to_dict().items():
            results = {}
            results["name"] = key
            results["y"] = round(value, 2)
            array_data.append(round(value, 2))
            if key == "Classified":
                results["sliced"] = True
                results["selected"] = True
            pie_chart_data.append(results)
        pie_chart_data.append(array_data)
        series = {"data": pie_chart_data, "name": "Classification Proportions", "colorByPoint": True}

        # drop unclassified from prop charts
        df = df[df["barcode_name"] != "unclassified"]
        listy = []
        for name, group in df.groupby("status"):
            listy.append({"name": name, "data": group["data"].values.tolist()})

        listy.append({"categories": categories, "pieChartData": series})

        return Response(listy, status=status.HTTP_200_OK)

    return Response("Not a barcoded run.", status=status.HTTP_204_NO_CONTENT)


@api_view(["GET"])
def flowcell_barcodes_list(request):
    """
    Return a list of barcodes for a given flowcell.
    Parameters
    ----------
    request: rest_framework.request.Request

    Returns
    -------

    """
    flowcell_id = request.GET.get("flowcellId", None)
    if not flowcell_id:
        return Response("No Flowcell ID specified.", status.HTTP_400_BAD_REQUEST)
    flowcell = get_object_or_404(Flowcell, pk=int(flowcell_id))
    # TODO this could become slow if we ended up with 10000s of barcodes
    barcodes = sorted(
        list(
            Barcode.objects.filter(run__in=flowcell.runs.all()).values_list(
                "name", flat=True
            ).distinct()
        )
    )

    if "Unblocked" not in barcodes:
        barcodes.remove("Sequenced")

    return Response(barcodes, status=status.HTTP_200_OK)


@api_view(["POST"])
def reactivate_flowcell(request, pk):
    """
    Reactivate a flowcell so it isn't ignored by celery run monitor
    Parameters
    ----------
    request: rest_framework.request.Request
    pk: int
        Primary key of the flowcell record to update

    Returns
    -------

    """

    flowcell = Flowcell.objects.get(pk=pk)
    if flowcell.active():
        return Response("Error - Flowcell already active", status=405)
    flowcell.last_activity_date = datetime.datetime.now(datetime.timezone.utc)
    flowcell.save()
    return Response("Flowcell last activity date successfully updated", status=200)


@api_view(["GET"])
def read_type_list(request):
    """
    :purpose: Returns a list of FastqReadType to the client
    :used_by: minotour client
    :author: Roberto Santos

    ChangeLog
    2018-07-09 Add documentation - Roberto

    :param request: (standard django request) without querystring parameter
    :return: (string) json format
    """
    queryset = FastqReadType.objects.all()
    serializer = FastqReadTypeSerializer(
        queryset, many=True, context={"request": request}
    )
    return Response(serializer.data)


@api_view(["GET"])
def read_type_detail(request, pk):
    """
    :purpose: Retrieve a FastqReadType instance
    :used_by: minotour client
    :author: Roberto Santos

    ChangeLog
    2018-07-09 Add documentation - Alex

    :param request: (standard django request) ???
    :param pk: (int) run primary key
    :return: ???
    """
    try:
        run = FastqReadType.objects.get(pk=pk)
    except FastqReadType.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    serializer = FastqReadTypeSerializer(run, context={"request": request})
    return Response(serializer.data)


@api_view(["GET"])
def events_type_list(request):
    """
    :purpose: Returns a list of all possible events to the client
    :used_by: minotour client
    :author: Roberto Santos

    ChangeLog
    2018-07-09 Add documentation - Rory

    :param request: (Django Request Object) No query parameters
    :return: (String) Json Format string of event types
    """
    queryset = MinionEventType.objects.all()
    serializer = MinIONEventTypeSerializer(
        queryset, many=True, context={"request": request}
    )
    return Response(serializer.data)


@api_view(["GET"])
def events_type_detail(request, pk):  # TODO consider removing
    """
    :purpose: Returns a single Minion Event to the client
    :used_by: minotour client
    :author: Roberto Santos

    ChangeLog
    2018-07-09 Add documentation - Rory

    :param request: (Django Request Object) No query parameters
    :param pk: (Integer) - Primary key for object lookup from database

    :return: (String) Json Format string of a single event types
    """

    try:
        event_ = MinionEventType.objects.get(pk=pk)
    except MinionEventType.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    serializer = MinIONEventTypeSerializer(event_, context={"request": request})
    return Response(serializer.data)


@api_view(["GET"])
def fastq_detail(request, pk):
    """

    :param request:
    :param pk: fastqid
    :return:
    """
    queryset = FastqFile.objects.filter(runid=pk)
    serializer = FastqFileSerializer(queryset, many=True, context={"request": request})
    return Response(serializer.data)


@api_view(["GET", "POST"])
def fastq_file(request, pk):
    """
    :purpose: returns an md5 checksum for a file as seen by minotour for a specific run id
    :used_by: minotour client gets data and checks observed files against it.
    :author: Matt Loose

    :param request: (standard django request) without querystring parameter
    :param pk: pk is the runid

    :return: (str) json format
    """
    """if request.method == "GET":
        try:
            queryset = FastqFile.objects.filter(runid=pk)
            serializer = FastqFileSerializer(queryset, many=True, context={'request': request})
            return Response(serializer.data)
        except FastqFile.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)
    """

    if request.method == "GET":

        queryset = FastqFile.objects.filter(runid=pk).filter(owner=request.user)

        serializer = FastqFileSerializer(
            queryset, many=True, context={"request": request}
        )

        return Response(serializer.data)

    elif request.method == "POST":

        run = (
            Run.objects.filter(runid=request.data["runid"])
                .filter(owner=request.user)
                .first()
        )

        obj, created = FastqFile.objects.get_or_create(
            runid=request.data["runid"], name=request.data["name"], owner=request.user,
        )
        obj.run = run
        obj.md5 = request.data["md5"]

        obj.save()

        serializer = FastqFileSerializer(obj, context={"request": request})

        return Response(serializer.data, status=status.HTTP_201_CREATED)

    return Response(status=status.HTTP_400_BAD_REQUEST)


@api_view(["GET", "POST"])
def run_list(request):
    """
    :purpose: Get returns a list of all runs in minotour owned by a specific user.
              Post will create a new run in minotour for a specific user.
    :used_by: minotour client gets and posts runs, minotour app gets runs only.
    :author: Roberto Santos

    ChangeLog
    2018-07-09 Add documentation - Matt

    :param request: (standard django request) without querystring parameter
    :return: (str) json format
    """

    if request.method == "GET":

        flowcell_list = []

        for flowcell in Flowcell.objects.all():
            if (
                    request.user == flowcell.owner
                    or request.user.has_perm("view_data", flowcell)
                    or request.user.has_perm("run_analysis", flowcell)
            ):
                flowcell_list.append(flowcell)

        print(flowcell_list)

        queryset = Run.objects.filter(flowcell__in=flowcell_list).filter(
            to_delete=False
        )
        serializer = RunSerializer(queryset, many=True, context={"request": request})

        return Response(serializer.data)

    elif request.method == "POST":

        serializer = RunSerializer(data=request.data, context={"request": request})

        if serializer.is_valid():
            run = serializer.save(owner=request.user)
            # Change the flowcell to active
            flowcell = run.flowcell

            flowcell.last_activity_date = datetime.datetime.now(datetime.timezone.utc)

            flowcell.save()
            # Create the rejected accepted barcodes
            # accepted_barcode = Barcode(run=run, name="S")

            # rejected_barcode = Barcode(run=run, name="U")

            # accepted_barcode.save()

            # rejected_barcode.save()

            return Response(serializer.data, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(["GET", "PUT"])
def run_detail(request, pk):
    """
    Retrieve, update or delete a run instance.
    """

    search_criteria = request.GET.get("search_criteria", "id")

    flowcell_list = []  # TODO move the next 3 lines to a function

    for flowcell in Flowcell.objects.all():
        if (
                request.user == flowcell.owner
                or request.user.has_perm("view_data", flowcell)
                or request.user.has_perm("run_analysis", flowcell)
        ):
            flowcell_list.append(flowcell)

    if search_criteria == "runid":

        run_list = Run.objects.filter(flowcell__in=flowcell_list).filter(runid=pk)

    elif search_criteria == "id":

        run_list = Run.objects.filter(flowcell__in=flowcell_list).filter(id=pk)

    else:

        run_list = Run.objects.none()

    if len(run_list) < 1:
        return Response(status=status.HTTP_404_NOT_FOUND)

    run = run_list[
        0
    ]  # TODO what if by any means we have more than one run if the same runid?

    if request.method == "GET":
        serializer = RunSerializer(run, context={"request": request})
        return Response(serializer.data)

    elif request.method == "PUT":
        serializer = RunSerializer(run, data=request.data, context={"request": request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(["GET", "POST"])
def minion_list(request):
    """
    Get returns a list of all minions for a given user. Post creates a new Minion for a given user.
    Parameters
    ----------
    request

    Returns
    -------

    """
    if request.method == "GET":
        queryset = Minion.objects.filter(owner=request.user)
        serializer = MinionSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    elif request.method == "POST":
        print(f"request data for minion POST is:")
        print(request.data)
        serializer = MinionSerializer(data=request.data, context={"request": request})
        if serializer.is_valid():
            serializer.save(owner=request.user)
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(["GET"])
def active_minion_list(request):
    """
    Get all active minions for a given user. Displays for the remote control page.
    Parameters
    ----------
    request: rest_framework.request.Request

    Returns
    -------
    A list of Active minions, with a dictionary of information about them.

    """
    active_minion_list = []

    extra_data = {}

    blank_stats = {
        "id": -1,
        "minion_id": -1,
        "run_id": -1,
        "sample_time": datetime.datetime.now(),
        "event_yield": "Unknown",
        "asic_temp": "Unknown",
        "heat_sink_temp": "Unknown",
        "voltage_value": "Unknown",
        "mean_ratio": "Unknown",
        "open_pore": "Unknown",
        "in_strand": "Unknown",
        "multiple": "Unknown",
        "unavailable": "Unknown",
        "unknown": "Unknown",
        "adapter": "Unknown",
        "pending_mux_change": "Unknown",
        "unclassified": "Unknown",
        "below": "Unknown",
        "unblocking": "Unknown",
        "above": "Unknown",
        "good_single": "Unknown",
        "saturated": "Unknown",
        "inrange": "Unknown",
        "strand": "Unknown",
        "no_pore": "Unknown",
        "zero": "Unknown",
        "pore": "Unknown",
        "minKNOW_read_count": "Unknown",
        "minKNOW_histogram_values": "[]",
        "minKNOW_histogram_bin_width": "Unknown",
        "created_date": datetime.datetime.now(),
        "actual_max_val": "Unknown",
    }

    data = {}

    for minion in Minion.objects.filter(owner=request.user):
        # Store extra data about actively sequencing minions here
        extra_data[minion.name] = {}

        minion_event_list = MinionEvent.objects.filter(minION=minion).order_by(
            "datetime"
        )

        extra_data[minion.name]["freespace"] = minion.space_available()
        extra_data[minion.name]["space_till_shutdown"] = minion.space_till_shutdown()
        extra_data[minion.name]["minknow_version"] = minion.minknow_version()

        if minion_event_list.count() > 0:

            last_minion_event = minion_event_list.last()
            five_minute_check = datetime.datetime.now(
                datetime.timezone.utc
            ) - datetime.timedelta(minutes=5)

            if (
                    last_minion_event.event.name != "unplugged"
                    and minion.currentdetails.last_modified > five_minute_check
            ):

                if minion.currentdetails.minKNOW_status in [
                    "ACQUISITION_RUNNING",
                    "ACQUISITION_PROCESSING",
                    "ACQUISITION_FINISHING",
                ]:
                    # Get the run start time and other deets
                    extra_data[minion.name]["start_time"] = str(minion.start_time())
                    extra_data[minion.name][
                        "wells_per_channel"
                    ] = minion.wells_per_channel()
                    extra_data[minion.name][
                        "read_length_type"
                    ] = minion.read_length_type()
                    extra_data[minion.name]["actual_max_val"] = minion.actual_max_val()
                    extra_data[minion.name]["target_temp"] = minion.target_temp()
                    extra_data[minion.name]["flowcell_type"] = minion.flowcell_type()
                    extra_data[minion.name][
                        "experiment_name"
                    ] = minion.experiment_name()

                    mrs = MinionRunStats.objects.filter(minion=minion).last()

                    if mrs:
                        mrs_serialiser = MinionRunStatsSerializer(
                            mrs, context={"request": request}
                        )

                        data[minion.name] = mrs_serialiser.data

                active_minion_list.append(minion)

    serializer = MinionSerializer(
        active_minion_list, many=True, context={"request": request}
    )

    return_data = list(serializer.data)

    for active_minion in return_data:
        active_minion.update(extra_data[active_minion["name"]])
        if active_minion.get("name", 0) in data:
            print("Hello")
            pprint(extra_data)
            active_minion.update(data.get(active_minion["name"], blank_stats))
            active_minion.update(extra_data[active_minion["name"]])

    pprint(serializer.data)
    return Response(serializer.data)


@api_view(["GET"])
def computers_list(request):
    """
    Get the computer information for remote control connected computers
    Parameters
    ----------
    request: rest_framework:request.Request
        Request body.
    Returns
    -------

    """


@api_view(["GET", "POST"])
def minion_messages_list(request, pk):
    """
    :purpose: Get returns a list of all messages sent to minotour from minKNOW for a specific MinION.
              Post creates a new message from a minION.
    :used_by: minotour app uses this endpoint to show messages from the MinION.
    :author: Roberto Santos

    ChangeLog
    2018-07-09 Add documentation

    :param request: (standard django request) without querystring parameter
    :param pk: (int) Minion id
    :return: (str) json format
    """

    if request.method == "GET":

        queryset = MinionMessage.objects.filter(minION=pk)

        serializer = MinionMessageSerializer(
            queryset, many=True, context={"request": request}
        )

        for message in serializer.data:
            if message["full_text"] == "":
                message["message"] = message["full_text"]

        return Response(serializer.data)

    elif request.method == "POST":
        if len(request.data["message"]) > 256:
            request.data["full_text"] = request.data["message"]

            request.data["message"] = request.data["message"][:256]

        serializer = MinionMessageSerializer(
            data=request.data, context={"request": request}
        )
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


# @api_view(["GET"],)
# def sinceminion_messages_list(request, pk, starttime, endtime):
#
#     correctedstart = parser.parse(starttime) - datetime.timedelta(minutes=180)
#
#     correctedend = parser.parse(endtime) + datetime.timedelta(minutes=180)
#
#     queryset = (
#         MinionMessage.objects.filter(minION=pk)
#         .filter(
#             minKNOW_message_timestamp__gte=correctedstart.isoformat().replace(
#                 "+00:00", "Z"
#             )
#         )
#         .filter(
#             minKNOW_message_timestamp__lte=correctedend.isoformat().replace(
#                 "+00:00", "Z"
#             )
#         )
#     )
#
#     serializer = MinionMessageSerializer(
#         queryset, many=True, context={"request": request}
#     )
#
#     return Response(serializer.data)


@api_view(["GET"], )
def recentminion_messages_list(request, pk):
    queryset = MinionMessage.objects.filter(minION=pk).filter(
        minKNOW_message_timestamp__gte=timezone.now() - datetime.timedelta(hours=24)
    )

    serializer = MinionMessageSerializer(
        queryset, many=True, context={"request": request}
    )

    return Response(serializer.data)


@api_view(["GET"])
def minknow_message_html(request, pk):
    """
    :purpose: Return messages generated by MinKnow for a particular date and time range.
    :used_by: Flowcell summary page at the web app.
    :author: Roberto Santos

    ChangeLog
    2018-07-09 Add documentation - Roberto

    :param request: (standard django request)
    :param request.start_time: (datetime) Start of the date and time range
    :param request.end_time: (datetime) End of the date and time range # TODO include format example
    :param pk: (int) Flowcell id
    :return: (str) Json format
    """
    flowcell = Flowcell.objects.get(pk=pk)
    runs = Run.objects.filter(flowcell=flowcell)

    messages = MinionMessage.objects.filter(run__in=runs).order_by("-timestamp")

    return render(request, "reads/minknow_messages.html", {"message_list": messages})


@api_view(["GET", "POST"], )
def minION_control_list(request, pk):
    """
    TODO describe function
    """
    if request.method == "GET":
        queryset = MinionControl.objects.filter(minION=pk).filter(complete=False)
        serializer = MinIONControlSerializer(
            queryset, many=True, context={"request": request}
        )
        return Response(serializer.data)

    elif request.method == "POST":
        serializer = MinIONControlSerializer(
            data=request.data, context={"request": request}
        )

        if serializer.is_valid():
            serializer.save(owner=request.user)
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(["GET", "POST"], )
def minION_control_update(request, pk, checkid):
    """
    TODO describe function
    """
    try:
        event_ = MinionControl.objects.get(id=checkid)

        if event_.owner != request.user:
            return Response(status=status.HTTP_401_UNAUTHORIZED)

    except MinionControl.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == "GET":
        serializer = MinIONControlSerializer(event_, context={"request": request})
        return Response(serializer.data)

    if request.method == "POST":
        event_ = MinionControl.objects.get(id=checkid)
        event_.complete = True
        event_.save(update_fields=["complete"])
        return Response(status=status.HTTP_204_NO_CONTENT)


@api_view(
    ["GET", ]
)
def minION_currentrun_list(request, pk):
    """
    TODO describe function
    """
    try:
        minion = Minion.objects.get(pk=pk)
    except Minion.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == "GET":
        serializer = MinionSerializer(minion, context={"request": request})
        return Response(serializer.data)


@api_view(["GET", "POST"])
def minION_events_list(request, pk):
    """
    List of all events associated with a minION or add a new event.
    """
    if request.method == "GET":
        queryset = MinionEvent.objects.filter(minION=pk)
        serializer = MinIONEventSerializer(
            queryset, many=True, context={"request": request}
        )
        return Response(serializer.data)

    elif request.method == "POST":
        serializer = MinIONEventSerializer(
            data=request.data, context={"request": request}
        )
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(["GET", "POST"])
def minION_scripts_list(request, pk):
    """
    List of all scripts associated with a minION.
    """
    if request.method == "GET":
        queryset = MinionScripts.objects.filter(minION=pk)
        serializer = MinIONScriptsSerializer(
            queryset, many=True, context={"request": request}
        )
        return Response(serializer.data)

    elif request.method == "POST":
        serializer = MinIONScriptsSerializer(
            data=request.data, partial=True, context={"request": request}
        )
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


# TODO THIS IS A CLIENT ENDPOINT
@api_view(["GET", "PUT", "POST", "DELETE"])
def list_minion_info(request, pk):
    """
    Update, Create, GET and Delete  minion info entries. Sent by the client, contains minKNOW and minION information.`
    Parameters
    ----------
    request: rest_framework.request.Request
        Request object from the client.
    pk: int
        primary key of the minion entry in the database.

    Returns
    -------
    Response: rest_framework.response.Response
        The response of the endpoint. May be an error, otherwise returns a copy of the minIONStatus.

    """

    if request.method == "POST":
        serializer = MinionInfoSerializer(
            data=request.data, context={"request": request}
        )
        if serializer.is_valid(raise_exception=True):
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        print(serializer.errors)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    try:
        minion_info = MinionInfo.objects.get(minion=pk)

    except MinionInfo.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == "GET":
        serializer = MinionInfoSerializer(minion_info, context={"request": request})
        return Response(serializer.data)

    elif request.method == "PUT":
        # TODO WE HAVE A BUG HERE AGAIN?
        serializer = MinionInfoSerializer(
            minion_info, data=request.data, context={"request": request}
        )
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    elif request.method == "DELETE":
        minion_info.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


@api_view(["GET", "PUT", "POST", "DELETE"])
def minion_run_stats_list(request, pk):
    # TODO CLIENT ONLY
    """
    Create, Update, Get or Delete minion_run_stats entries in the database. Stores run INFO that is live fed to us by MinKNOW/minFQ.
    Parameters
    ----------
    request: rest_framework.request.Request
        The request object of AJAX request form the client.
    pk: str
        The url instance that links to the run in the Database. Provided by client.

    Returns
    -------

    """
    if request.method == "POST":

        serializer = MinionRunStatsSerializer(
            data=request.data, context={"request": request}
        )
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    try:
        crazy_minion_run_stats = MinionRunStats.objects.filter(run_id=pk)

    except MinionRunStats.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == "GET":

        serializer = MinionRunStatsSerializer(
            crazy_minion_run_stats, many=True, context={"request": request}
        )

        return Response(serializer.data)

    elif request.method == "PUT":
        serializer = MinionRunStatsSerializer(
            crazy_minion_run_stats, data=request.data, context={"request": request}
        )
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    elif request.method == "DELETE":
        crazy_minion_run_stats.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


@api_view(["GET"])
def minion_run_stats_latest(request, pk, checkid):
    """
    TODO Web
    Get the latest minion RUN stats, to be displayed in the Live event data page.
    Parameters
    ----------
    request: rest_framework.request.Request
        Ajax request body.
    pk: int
        Primary key (ID) of the run that these stats pertain to.
    checkid: int
        Id of last entry that we checked.

    Returns
    -------
    list
        A list of dictionaries of the Live Event data stats.
    """

    try:
        crazy_minion_run_stats = MinionRunStats.objects.filter(
            run_id=pk, id__gt=checkid
        )[:1000]

    except MinionRunStats.DoesNotExist:

        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == "GET":
        serializer = MinionRunStatsSerializer(
            crazy_minion_run_stats, many=True, context={"request": request}
        )

        return Response(serializer.data)


@api_view(["GET"])
def remote_control_run_stats(request):
    """
    Get the latest run stats for a run, for a given user.
    Parameters
    ----------
    request: rest_framework.request.Request
        The Ajax request body.
    Returns
    -------
    data: list
        List of dictionaries. Key minion device ID, value run stats.
    """

    data = {}

    minions = request.GET.get("activeMinions", [])

    if minions:
        minions = json.loads(minions)

        for minion in minions:

            queryset = MinionRunStats.objects.filter(
                minion__owner=request.user, minion__name=minion
            ).values()

            if queryset:
                last_runs_stat = queryset[::-1][0]

                data[minion] = last_runs_stat

        return Response(data, status=200)

    return Response(status=400)


@api_view(["GET", "PUT", "POST", "DELETE"])
def list_minion_run_status(request, pk):
    """
    TODO Client and Web
    Get, Create, Update or Delete MinionRunINfo records from the Minotour Database.
    Parameters
    ----------
    request: rest_framework.request.Request
        Request object from AJAX request, both client and web.
    pk: int
        The Pk of the Run that these Metadata are pertaining to.

    Returns
    -------

    """
    print(request.data)
    print("minion tun info")
    if request.method == "POST":
        serializer = MinionRunInfoSerializer(
            data=request.data, context={"request": request}
        )
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    try:
        minion_run_info = MinionRunInfo.objects.filter(run_id=pk)

    except MinionRunInfo.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == "GET":
        serializer = MinionRunInfoSerializer(
            minion_run_info, many=True, context={"request": request}
        )
        return Response(serializer.data)

    elif request.method == "PUT":
        serializer = MinionRunInfoSerializer(
            minion_run_info, data=request.data, context={"request": request}
        )
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    elif request.method == "DELETE":
        minion_run_info.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


@api_view(["GET", "PUT", "DELETE"])
def minION_scripts_detail(request, pk, nk):
    """
    TODO describe function
    """
    try:
        script = MinionScripts.objects.get(pk=nk)
    except MinionScripts.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == "GET":
        serializer = MinIONScriptsSerializer(script, context={"request": request})
        return Response(serializer.data)

    elif request.method == "PUT":
        serializer = MinIONScriptsSerializer(
            script, data=request.data, context={"request": request}
        )
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    elif request.method == "DELETE":
        script.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


# @api_view(["GET", "POST"])
# def read_list(request, pk):
#     if request.method == "GET":
#
#         querysets = FastqRead.objects.filter(run_id=pk)
#
#         serializer = FastqReadGetSerializer(
#             querysets, many=True, context={"request": request}
#         )
#
#         return Response(serializer.data)
#
#     elif request.method == "POST":
#
#         # serializer = FastqReadSerializer(
#         #     data=request.data, many=True, context={"request": request}
#         # )
#         #
#         # if serializer.is_valid():
#         #     serializer.save()
#         #     return Response(serializer.data, status=status.HTTP_201_CREATED)
#
#         reads = request.data
#         reads_list = [FastqRead(**read) for read in reads]
#         try:
#             FastqRead.objects.bulk_create(read_list)
#             return Response("Updated reads", status=status.HTTP_201_CREATED)
#         except Exception as e:
#             return Response(e, status=status.HTTP_400_BAD_REQUEST)


@api_view(["GET", "POST"])
def minION_detail(request, pk):
    # """
    # TODO describe function
    # """
    print(request.user)
    print(type(request.user))
    print(request.GET.get("search_criteria", "LOLNO"))
    print(pk)

    if request.method == "GET":

        search_criteria = request.GET.get("search_criteria", "id")

        if search_criteria == "id":

            minion_list = Minion.objects.filter(owner=request.user, id=pk)

        elif search_criteria == "name":

            minion_list = Minion.objects.filter(owner=request.user, name=pk)

        else:

            minion_list = Minion.objects.none()

        print(minion_list)

        if len(minion_list) == 0:
            return Response(status=status.HTTP_404_NOT_FOUND)

        minion = minion_list[0]

        serializer = MinionSerializer(minion, context={"request": request})

        return Response(serializer.data)

    elif request.method == "POST":

        minion = Minion.objects.get(pk=pk)

        serializer = MinionSerializer(
            minion, data=request.data, partial=True, context={"request": request}
        )

        if serializer.is_valid():
            serializer.save()

        else:
            print("PROBLEM")

        return Response(serializer.data)


@api_view(["GET"])
def readname_list(request, pk):
    """
    TODO describe function
    CUrrent API behaviour asssumes starting from page1.
    Request for page 0 returns last page.
    This could be a source of confusion and we should resolve.
    """
    if request.method == "GET":
        queryset = FastqRead.objects.filter(fastqfile_id=pk).order_by("id")

        paginator = Paginator(queryset, settings.PAGINATION_PAGE_SIZE)

        page = request.GET.get("page")

        try:
            result = paginator.page(page)
        except PageNotAnInteger:
            result = paginator.page(1)
        except EmptyPage:
            result = paginator.page(paginator.num_pages)

        result2 = set()
        for key in result:
            result2.add(key.read_id)

        result3 = {}
        result3["number_pages"] = paginator.num_pages
        result3["data"] = list(result2)
        return HttpResponse(json.dumps(result3), content_type="application/json")


@api_view(["GET"])
def run_summary_barcode(request, pk):
    """
    Return a list with summaries for each read type of a given run.
    """
    queryset = RunSummaryBarcode.objects.filter(run_id__owner=request.user).filter(
        run_id=pk
    )

    serializer = RunSummaryBarcodeSerializer(
        queryset, many=True, context={"request": request}
    )

    return Response(serializer.data)


@api_view(["GET"])
def run_summary_barcode_by_minute(request, pk):
    """
    Return a list with summaries for a particular run grouped by minute.
    """
    queryset = RunStatisticBarcode.objects.filter(run_id__owner=request.user).filter(
        run_id=pk
    )

    serializer = RunStatisticBarcodeSerializer(
        queryset, many=True, context={"request": request}
    )

    return Response(serializer.data)


@api_view(["GET"])
def run_channel_summary(request, pk):
    """
    Return a list with channel info for a particular run.
    """
    queryset = ChannelSummary.objects.filter(run_id__owner=request.user).filter(
        run_id=pk
    )

    serializer = ChannelSummarySerializer(
        queryset, many=True, context={"request": request}
    )

    return Response(serializer.data)


@api_view(["GET", "POST"])
def barcode_list(request, pk):
    """
    Return a list of barcodes of a particular run.
    """

    if request.method == "GET":
        queryset = Barcode.objects.filter(run__owner=request.user).filter(run__id=pk)

        serializer = BarcodeSerializer(
            queryset, many=True, context={"request": request}
        )

        return Response(serializer.data)

    elif request.method == "POST":

        barcodename = request.data["name"]

        minionrun = Run.objects.get(pk=pk)

        barcode, created2 = Barcode.objects.get_or_create(
            run=minionrun, name=barcodename
        )

        serializer = BarcodeSerializer(barcode, context={"request": request})

        return Response(serializer.data, status=status.HTTP_201_CREATED)


@api_view(["GET"])
def barcode_detail(request, pk):
    """
    Return the details of a particular barcode.
    """
    queryset = Barcode.objects.get(pk=pk)

    serializer = BarcodeSerializer(queryset, many=False, context={"request": request})

    return Response(serializer.data)


@api_view(["GET", "POST"])
def flowcell_list(request):
    """
    API endpoint for either GETting a list of all the flowcells for a given user (The request maker)
    Or for POSTing a new flowcells data to create a database entry
    :param request: Django rest framework request object
    :return:  A list of each flowcell details in dictionary form
    """
    if request.method == "GET":

        flowcells = []

        # queryset = Flowcell.objects.filter(owner=request.user)
        queryset = Flowcell.objects.all()

        for record in queryset:

            if request.user == record.owner:
                owner = True
                permission = "Owner"
            elif request.user.has_perm("view_data", record):
                owner = False
                permission = "View data"
            elif request.user.has_perm("run_analysis", record):
                owner = False
                permission = "Run analysis"
            else:
                logger.error(
                    f"Flowcell {record.name} retrieved for user {request.user} that should not have been."
                )
                continue

            obj = {
                "id": record.id,
                "name": record.name,
                "size": record.size,
                "start_time": record.start_time,
                "number_reads": record.number_reads,
                "number_reads_processed": record.number_reads_processed,
                "number_runs": record.number_runs,
                "number_barcodes": record.number_barcodes,
                "total_read_length": record.total_read_length,
                "average_read_length": record.average_read_length,
                "is_active": record.active(),
                "sample_name": record.sample_name,
                "has_fastq": record.has_fastq,
                "owner": owner,
                "permission": permission,
            }

            flowcells.append(obj)

        return Response({"data": flowcells})
    # If the method is not GET, it must be post
    else:
        # Serialise the data inside the request
        serializer = FlowcellSerializer(data=request.data, context={"request": request})
        # If the serialiser is valid
        if serializer.is_valid():
            # Save the data, creating the new flowcell
            serializer.save(owner=request.user)
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        # If the serialiser data is not valid, return bad request
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(["GET", "POST"])
def flowcell_detail(request, pk):
    """
    Get the details for a single flowcell
    :param request: The django rest framework request object
    :type request: rest_framework.request.Request
    :param pk: The primary key of the flowcell database row
    :type pk: int
    :return:
    """
    if request.method == "GET":

        search_criteria = request.GET.get("search_criteria", "id")

        if search_criteria == "id":
            pk = int(pk)

            flowcell_list = Flowcell.objects.filter(owner=request.user, id=pk)

        elif search_criteria == "name":

            flowcell_list = Flowcell.objects.filter(owner=request.user, name=pk)
        else:

            flowcell_list = Flowcell.objects.none()

        if not flowcell_list:

            # flowcell_list = return_shared_flowcells(pk, request)
            flowcell_list = []

            all_flowcell_list = Flowcell.objects.exclude(owner=request.user)

            for f in all_flowcell_list:

                if request.user.has_perm("view_data", f) or request.user.has_perm(
                        "run_analysis", f
                ):

                    if search_criteria == "id":
                        pk = int(pk)
                        if f.id == pk:
                            flowcell_list.append(f)
                    elif search_criteria == "name":
                        if f.name == pk:
                            flowcell_list.append(f)

        # TODO look at the logic here, a user with two identical flowcells would crash this

        if len(flowcell_list) != 1:
            return Response({"data": {}}, status=404)

        # get the matching flowcell
        flowcell = flowcell_list[0]
        # deserialise the flowcell data
        serializer = FlowcellSerializer(flowcell, context={"request": request})

        data = serializer.data
        # Add in the last read start time, to filter the MinKnow messages inside the run time
        # Enumerate through the data["runs"] value, an array containing a dict element, one for each run
        for i, array in enumerate(data["runs"]):
            # Get the last read time from the run summary by doing FK relationship lookups
            run = flowcell.runs.get(runid=data["runs"][i]["runid"])

            if hasattr(run, "summary"):
                last_read_time = run.summary.last_read_start_time
            else:
                last_read_time = datetime.datetime.now(datetime.timezone.utc)
            # Add it to the dictionary that contains information on the run
            data["runs"][i]["last_read_time"] = last_read_time

        return_dict = {"data": data}

        return Response(return_dict)

    elif request.method == "POST":

        serializer = FlowcellSerializer(data=request.data, context={"request": request})

        if serializer.is_valid():
            serializer.save()

            return Response(serializer.data, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(["GET"])
def flowcell_summary_barcode(request, pk):
    q = request.GET.get("barcode_name", "All reads")

    flowcell = Flowcell.objects.get(pk=pk)

    qs = FlowcellSummaryBarcode.objects.filter(flowcell=flowcell).filter(
        barcode_name=q
    )  # .filter(run_id__owner=request.user)

    serializer = FlowcellSummaryBarcodeSerializer(qs, many=True)

    return Response(serializer.data)


def process_summary_data(df):
    """
    This function expects a total_length, read_count and quality_sum in the dataframe.
    :param df: dataframe to add percentage, avg length and avq quality
    :return: dataframe
    """

    df["percentage"] = np.where(
        df["barcode_name"] == "All reads",
        df["total_length"]
        / df[df["barcode_name"].eq("All reads")]["total_length"].sum()
        * 100,
        df["total_length"]
        / df[df["barcode_name"].ne("All reads")]["total_length"].sum()
        * 100,
    )

    df["avg_length"] = df["total_length"] / df["read_count"]

    df["avg_quality"] = df["quality_sum"] / df["read_count"]

    return df


@login_required
def flowcell_basecalled_summary_html(request, pk):
    """
    Returns a html table to the top of the chancalc page
    :param request: The Get request body
    :type request: rest_framework.request.Request
    :param pk: Primary key of the flowcell
    :type pk: int
    :return: application/html table
    """

    flowcell = Flowcell.objects.get(pk=pk)
    # Get the flowcell summary barcodes for this flowcell, generated during the task
    qs = FlowcellSummaryBarcode.objects.filter(flowcell=flowcell).exclude(
        barcode_name="No barcode"
    )
    # Create a dataframe of the fetched data
    df = pd.DataFrame.from_records(
        qs.values(
            "barcode_name",
            "rejection_status",
            "max_length",
            "min_length",
            "quality_sum",
            "read_count",
            "read_type_name",
            "status",
            "total_length",
        )
    )
    # TODO do we not already keep these - do this on creation of flowcell summary barcode?!?!
    # calculate the percentage of reads in a barcode
    df = process_summary_data(df)

    df["read_type"] = np.where(df["status"] == "True", "Pass", "Fail")
    # remove status column
    df = df.drop("status", axis=1)
    # Create a record oriented dictionary from the dataframe
    dictdf = df.to_dict("records")
    # aggregations to be calculated
    agg = {
        "max_length": "max",
        "min_length": "min",
        "quality_sum": "sum",
        "read_count": "sum",
        "total_length": "sum",
    }
    # Group by All reads and pass fail read type, perform the above aggregations
    # and creating the columns for the results
    df2 = (
        df[df["barcode_name"].eq("All reads")]
            .groupby(["barcode_name", "read_type_name", "rejection_status"])
            .agg(agg)
    )

    df2 = df2.reset_index()

    df2 = process_summary_data(df2)

    df2["read_type"] = "all"
    # Append results to the dictionary
    dictdf.extend(df2.to_dict("records"))
    # Same aggregations on other barcodes
    df2 = (
        df[df["barcode_name"].ne("All reads")]
            .groupby(["barcode_name", "read_type_name", "rejection_status"])
            .agg(agg)
    )

    if df2.shape[0] > 0:
        df2 = df2.reset_index()

        df2 = process_summary_data(df2)

        df2["read_type"] = "all"

        dictdf.extend(df2.to_dict("records"))

    return render(request, "reads/flowcell_summary.html", {"qs": dictdf})


@api_view(["GET"])
def flowcell_statistics(request, pk):
    """
    Return a prepared set of summaries for reads quality over time grouped by minute. Used in the
    rest of the charts on the basecalled data tab - except pore charts
    :param request: The Get request body
    :type request: rest_framework.request.Request
    :param pk: Primary key of the flowcell the read data is attached to
    :type pk: int
    :return: (HttpResponse) json of all sorts of data
    """
    # TODO comment and REWRITE OOFT
    barcode_name = request.GET.get("barcodeName", "All reads")

    barcodes_list = [barcode_name]

    if barcode_name != "All reads":
        barcodes_list.append("All reads")

    flowcell = Flowcell.objects.get(pk=pk)
    # Reverse lookup a queryset of all the runs in this flowcell
    run_list = flowcell.runs.all()

    """
    queryset = (
        FlowcellStatisticBarcode.objects.filter(flowcell=flowcell, iteration_count=job_master_iteration)
            .filter(
            Q(barcode_name__in=barcodes_list) | Q(rejection_status__in=barcodes_list)
        )
            .order_by("sample_time")
    )
    """
    queryset = (
        FlowcellStatisticBarcode.objects.filter(flowcell=flowcell)
            .filter(
            Q(barcode_name__in=barcodes_list) | Q(rejection_status__in=barcodes_list)
        )
            .order_by("sample_time")
    )

    df = pd.DataFrame(
        list(
            queryset.values(
                "read_type_name",
                "sample_time",
                "barcode_name",
                "status",
                "read_count",
                "max_length",
                "total_length",
                "quality_sum",
                "rejection_status",
            )
        )
    )

    df["read_type"] = np.where(df["status"] == "True", "Pass", "Fail")

    agg = {
        "max_length": "max",
        "quality_sum": "sum",
        "read_count": "sum",
        "total_length": "sum",
    }

    df2 = df.groupby(
        ["barcode_name", "read_type_name", "sample_time", "rejection_status"]
    ).agg(agg)

    df = df.drop("status", axis=1)

    df2["read_type"] = "All"

    df2 = df2.reset_index()

    df3 = pd.concat([df, df2], ignore_index=True, sort=True)

    # TODO np.cumsum()
    df3["cumulative_read_count"] = df3.groupby(
        ["barcode_name", "read_type_name", "read_type", "rejection_status"]
    )["read_count"].apply(lambda x: x.cumsum())

    df3["cumulative_bases"] = df3.groupby(
        ["barcode_name", "read_type_name", "read_type", "rejection_status"]
    )["total_length"].apply(lambda x: x.cumsum())

    df3["key"] = (
            df3["barcode_name"].astype("str")
            + " - "
            + df3["read_type_name"].astype("str")
            + " - "
            + df3["read_type"].astype("str")
            + " - "
            + df3["rejection_status"].astype("str")
    )

    df3["average_quality"] = (
        df3["quality_sum"].div(df3["read_count"]).astype("float").round(decimals=2)
    )

    df3["average_quality"] = df3["average_quality"].astype("float")

    df3["average_length"] = df3["total_length"].div(df3["read_count"]).round(decimals=0)

    df3["sequence_rate"] = df3["total_length"].div(60).round(decimals=0)

    df3["corrected_time"] = df3["sample_time"].astype(np.int64) // 10 ** 6

    if barcode_name != "All reads":
        df3 = df3.drop(
            df3.index[(df3.barcode_name == "All reads") & (df3.read_type != "All")]
        )

    data_keys = df3["key"].unique().tolist()

    result_dict = {
        k: df3[df3["key"].eq(k)][
            [
                "corrected_time",
                "average_quality",
                "average_length",
                "cumulative_bases",
                "cumulative_read_count",
                "max_length",
                "sequence_rate",
            ]
        ].values.tolist()
        for k in df3.key.unique()
    }
    run_data = []

    for run in run_list:
        run_dict = {
            "id": run.id,
            "name": run.name,
            "runid": run.runid,
            "start_time": run.start_time,
            "last_read": run.last_read(),
        }

        run_data.append(run_dict)

    res = {
        "runs": run_data,
        "data_keys": data_keys,
        "data": result_dict,
        "date_created": datetime.datetime.now(),
    }

    return HttpResponse(
        json.dumps(res, cls=DjangoJSONEncoder), content_type="application/json"
    )


# @api_view(["GET"])
# def flowcell_speed(request, pk):
#     window = 5
#     flowcell = Flowcell.objects.get(pk=pk)
#
#     queryset = (
#         FlowcellStatisticBarcode.objects.filter(flowcell=flowcell)
#         .filter(barcode_name="All Reads")
#         .order_by("sample_time")
#     )
#
#     df = pd.DataFrame(
#         list(
#             queryset.values("sample_time", "status", "channel_presence", "total_length")
#         )
#     )
#     if len(df.status.unique()) == 2:
#         # We have pass and fail reads in the dataset
#         # Convert channel presence strings to numpy array for later manipulation.
#         df["channel_presence"] = df["channel_presence"].apply(
#             lambda x: np.asarray(list(x), dtype=int)
#         )
#         # Convert the dataframe to stack channel_presence and total_length
#         df = (
#             df.groupby("sample_time")[["channel_presence", "total_length"]]
#             .apply(lambda df: df.reset_index(drop=True))
#             .unstack()
#         )
#         # Generate the total length across both pass and fail
#         df["sum_total_length"] = df["total_length"][0] + df["total_length"][1]
#         # Add the two channel presence arrays to obtain 0 values and values greater than 1 where a channel has been used
#         df["sum_channel_presence"] = (
#             df["channel_presence"][0] + df["channel_presence"][1]
#         )
#         # Calculate our channel count by counting the number of nonzero values in the array.
#         # This gets around the fact that a channel might produce both a pass and fail read in the same one minute window
#         df["chan_count"] = df["sum_channel_presence"].apply(
#             lambda x: np.count_nonzero(x)
#         )
#         # calculate the mean channel count over a 10 minute window.
#         df["mean_chan_count"] = df["chan_count"].rolling(window=window).mean()
#         # calculate the mean read length seen in a rolling 10 minute window.
#         df["mean_total_length"] = df["sum_total_length"].rolling(window=window).mean()
#         # calculate the mean speed over those rolling windows in bases per second
#         df["mean_speed"] = (
#             df["mean_total_length"]
#             .div(df["mean_chan_count"])
#             .div(600)
#             .round(decimals=0)
#         )
#     elif len(df.status.unique()) == 1:
#         # We only have either pass or fail reads in the dataset
#         df["channel_presence"] = df["channel_presence"].apply(
#             lambda x: np.asarray(list(x), dtype=int)
#         )
#         df = df.groupby("sample_time")[["channel_presence", "total_length"]].apply(
#             lambda df: df.reset_index(drop=True)
#         )
#         df["chan_count"] = df["channel_presence"].apply(lambda x: np.count_nonzero(x))
#         df["mean_chan_count"] = df["chan_count"].rolling(window=window).mean()
#         df["mean_total_length"] = df["total_length"].rolling(window=window).mean()
#         df["mean_speed"] = (
#             df["mean_total_length"].div(df["mean_chan_count"]).round(decimals=0)
#         )
#         df.reset_index(level=1, drop=True, inplace=True)
#
#     return Response(df["mean_speed"].to_json(orient="columns"))


@api_view(["GET"])
def flowcell_histogram_summary(request, pk):
    """
    Return the data for histograms on the basecalled data page
    :param request: Get request body
    :type request: rest_framework.request.Request
    :param pk: The primary key of the flowcell that the basecalled data metrics are attached to
    :type pk: int
    :return: (rest_framework.response.Response) containing a
    """
    # Get the barcode we are doing this for, default to all reads
    barcode_name = request.GET.get("barcodeName", "All reads").replace("_", " ")
    # Get the information for the flowcell the data is attached to
    flowcell = Flowcell.objects.get(pk=pk)

    initialise = request.GET.get("initialise", False)

    total_reads_length = flowcell.total_read_length
    total_reads_count = flowcell.number_reads

    max_bin_index = (
        FlowcellHistogramSummary.objects.filter(flowcell=flowcell)
            .filter(Q(barcode_name=barcode_name) | Q(rejection_status=barcode_name))
            .aggregate(Max("bin_index"))
    )
    # create the categories for the x asis, they are the bins * the bin width (default 1000)
    categories = list(
        range(
            1 * FlowcellHistogramSummary.BIN_WIDTH,
            (max_bin_index["bin_index__max"] + 1) * FlowcellHistogramSummary.BIN_WIDTH,
            FlowcellHistogramSummary.BIN_WIDTH,
        )
    )

    # We are initialising the charts we only want the categories
    if initialise:
        return Response(categories, status=status.HTTP_200_OK)

    # A list of tuples, distinct on barcode name, pass/fail, template
    key_list = (
        FlowcellHistogramSummary.objects.filter(flowcell=flowcell)
            .filter(Q(barcode_name=barcode_name) | Q(rejection_status=barcode_name))
            .values_list("barcode_name", "read_type_name", "rejection_status", "status")
            .distinct()
    )

    chart_data = defaultdict(list)
    # unpack the tuple and create data using it, for each tuple
    for (barcode_name, read_type_name, rejection_status, is_pass) in key_list:
        # query flowcell Histogram summary objects for all of the combinations under this barcode
        histogram_queryset = (
            FlowcellHistogramSummary.objects.filter(flowcell=flowcell)
                .filter(Q(barcode_name=barcode_name) | Q(rejection_status=barcode_name))
                .filter(read_type_name=read_type_name)
                .filter(status=is_pass)
                .filter(rejection_status=rejection_status)
                .order_by("bin_index")
        )
        # make a list of 0s one element for each bin
        result_read_count_sum = np.zeros(
            max_bin_index["bin_index__max"] + 1, dtype=int
        ).tolist()
        # duplicate it into a new copy
        result_read_length_sum = result_read_count_sum.copy()

        result_collect_read_count_sum = [0] * (max_bin_index["bin_index__max"] + 1)

        result_collect_read_length_sum = [0] * (max_bin_index["bin_index__max"] + 1)
        # if there is more than one histogram summary, there should be one for each bin
        if histogram_queryset.count() > 0:
            if is_pass == "True":
                summary_is_pass = "Pass"
            else:
                summary_is_pass = "Fail"

            series_name = "{} - {} - {} - {}".format(
                barcode_name, read_type_name, rejection_status, summary_is_pass
            )

            seriescollect_read_count_sum = 0

            seriescollect_read_length_sum = 0

            for bin_summary in histogram_queryset:
                seriesbin_index = bin_summary.bin_index

                seriesread_count_sum = bin_summary.read_count

                seriesread_length_sum = bin_summary.read_length

                seriescollect_read_count_sum += bin_summary.read_count

                seriescollect_read_length_sum += bin_summary.read_length

                result_read_count_sum[seriesbin_index] = seriesread_count_sum

                result_read_length_sum[seriesbin_index] = seriesread_length_sum

                result_collect_read_count_sum[
                    seriesbin_index
                ] = seriescollect_read_count_sum

                result_collect_read_length_sum[
                    seriesbin_index
                ] = seriescollect_read_length_sum

            result_collect_read_count_sum = list(
                (
                        -pd.concat(
                            [
                                pd.Series([0]),
                                pd.Series(result_collect_read_count_sum).replace(
                                    to_replace=0, method="ffill"
                                ),
                            ]
                        )
                        + pd.Series(result_collect_read_count_sum)
                        .replace(to_replace=0, method="ffill")
                        .max()
                )
                / total_reads_count
                * 100
            )

            result_collect_read_length_sum = list(
                (
                        -pd.concat(
                            [
                                pd.Series([0]),
                                pd.Series(result_collect_read_length_sum).replace(
                                    to_replace=0, method="ffill"
                                ),
                            ]
                        )
                        + pd.Series(result_collect_read_length_sum)
                        .replace(to_replace=0, method="ffill")
                        .max()
                )
                / total_reads_length
                * 100
            )

            chart_data["read_count"].append(
                {"name": series_name, "data": result_read_count_sum}
            )
            chart_data["read_length"].append(
                {"name": series_name, "data": result_read_length_sum}
            )
            chart_data["collect_read_count"].append(
                {"name": series_name, "data": result_collect_read_count_sum, }
            )
            chart_data["collect_read_length"].append(
                {"name": series_name, "data": result_collect_read_length_sum, }
            )

    return Response({"data": chart_data, "categories": categories})


@api_view(["GET"])
def flowcell_channel_summary(request, pk):
    """
    This function generates the json data for the heatmap charts
    (Reads per channel and Bases per channel). The javascript
    file basecalledDataController.js calls this function.
    """

    flowcell = Flowcell.objects.get(pk=pk)

    qs = FlowcellChannelSummary.objects.filter(flowcell=flowcell)

    #
    # set all channels with value 0
    #

    result_mapped_to_flowcell = {}

    for i in range(1, flowcell.size + 1):
        coordinate = get_coords(i, flowcell.size)

        result_mapped_to_flowcell[i] = [coordinate[0], coordinate[1], 0, 0]

    #
    # fill the channels with correct value
    #

    for record in qs:
        coordinate = get_coords(record.channel, flowcell.size)

        result_mapped_to_flowcell[record.channel] = [
            coordinate[0],
            coordinate[1],
            record.read_count,
            record.read_length,
        ]

    return HttpResponse(
        json.dumps(list(result_mapped_to_flowcell.values()), cls=DjangoJSONEncoder),
        content_type="application/json",
    )


@api_view(["GET"])
def flowcell_run_summaries_html(request, pk):
    flowcell = Flowcell.objects.get(pk=pk)

    result = []

    for run in flowcell.runs.all():

        # Create an element with information from the run (originally from the fastq files)

        element = {
            "runid": None,
            "run_start_time": None,
            "first_read": None,
            "last_read": None,
            "minknow_computer_name": None,
            "minion_id": None,
            "asic_id": None,
            "sequencing_kit": None,
            "purpose": None,
            "minknow_version": None,
            "flowcell_type": None,
            "flowcell_id": None,
            "sample_name": None,
            "experiment_name": None,
        }

        minion_run_status_list = MinionRunInfo.objects.filter(run_id=run).order_by(
            "minKNOW_start_time"
        )

        if len(minion_run_status_list) > 0:
            minion_run_status = minion_run_status_list[0]

            element["runid"] = minion_run_status.run.runid
            element["run_start_time"] = minion_run_status.minKNOW_start_time
            element["minknow_computer_name"] = minion_run_status.minKNOW_computer
            element["minion_id"] = minion_run_status.minion.minION_name
            element["asic_id"] = minion_run_status.minKNOW_asic_id
            element["sequencing_kit"] = minion_run_status.sequencing_kit
            element["purpose"] = minion_run_status.minKNOW_exp_script_purpose
            element["minknow_version"] = minion_run_status.minKNOW_version
            element["flowcell_type"] = minion_run_status.flowcell_type
            element["flowcell_id"] = minion_run_status.minKNOW_flow_cell_id
            element["sample_name"] = minion_run_status.minKNOW_sample_name
            element["experiment_name"] = minion_run_status.experiment_id

            result.append(element)

    return render(request, "reads/flowcell_runs_summary.html", {"run_list": result})


@api_view(["GET"])
def flowcell_run_basecalled_summary_html(request, pk):
    flowcell = Flowcell.objects.get(pk=pk)

    result_basecalled_summary = []

    for run in flowcell.runs.all():

        if hasattr(run, "summary"):

            run_summary = {
                "runid": run.summary.run.runid,
                "read_count": run.summary.read_count,
                "total_read_length": run.summary.total_read_length,
                "max_read_length": run.summary.max_read_length,
                "min_read_length": run.summary.min_read_length,
                "avg_read_length": run.summary.avg_read_length,
                "first_read_start_time": run.summary.first_read_start_time,
                "last_read_start_time": run.summary.last_read_start_time,
            }

        else:

            print("RunSummary does not exist for run {}".format(run.id))

            run_summary_obj = return_temp_empty_summary(run)

            run_summary = {
                "runid": "Unavailable",
                "read_count": "Unavailable",
                "total_read_length": "Unavailable",
                "max_read_length": "Unavailable",
                "min_read_length": "Unavailable",
                "avg_read_length": "Unavailable",
                "first_read_start_time": "Unavailable",
                "last_read_start_time": "Unavailable",
            }

        result_basecalled_summary.append(run_summary)

    return render(
        request,
        "reads/flowcell_run_basecalled_summary.html",
        context={"result_basecalled_summary": result_basecalled_summary},
    )


@api_view(["GET"])
def flowcell_run_stats_latest(request, pk, checkid):
    """
    API endpoint for getting the live events chart data.
    :param request: DRF request body
    :param pk: The flowcell ID
    :param checkid: Check ID is defaulted to 0 at this time
    :return:
    """
    flowcell = Flowcell.objects.get(pk=pk)
    # Get all the Minion Run status objects for this flowcell
    minion_run_status_list = MinionRunInfo.objects.filter(
        run_id__flowcell=flowcell
    ).exclude(experiment_type="platform_qc")
    # If we have more than 0
    if minion_run_status_list.count() > 0:
        # Take the first entry to be the run status
        minion_run_status = minion_run_status_list[0]

    else:

        minion_run_status = None

    # Temporary work around to show all data.
    crazy_minion_run_stats = MinionRunStats.objects.filter(
        run_id__in=flowcell.runs.all().exclude(name="mux scan"), id__gt=checkid
    )

    result = []

    for minion_run_stats in crazy_minion_run_stats:
        element = {
            "id": minion_run_stats.id,
            # "minION": None,
            # "run_id": None,
            "sample_time": minion_run_stats.sample_time,
            "event_yield": minion_run_stats.event_yield,
            "asic_temp": round(minion_run_stats.asic_temp, 2),
            "heat_sink_temp": round(minion_run_stats.heat_sink_temp, 2),
            "voltage_value": minion_run_stats.voltage_value,
            # "mean_ratio": minion_run_stats.mean_ratio,
            "open_pore": minion_run_stats.open_pore,
            "in_strand": minion_run_stats.in_strand,
            "multiple": minion_run_stats.multiple,
            "unavailable": minion_run_stats.unavailable,
            "unknown": minion_run_stats.unknown,
            "adapter": minion_run_stats.adapter,
            "pending_mux_change": minion_run_stats.pending_mux_change,
            "unclassified": minion_run_stats.unclassified,
            "below": minion_run_stats.below,
            "unblocking": minion_run_stats.unblocking,
            "above": minion_run_stats.above,
            "good_single": minion_run_stats.good_single,
            "saturated": minion_run_stats.saturated,
            "inrange": minion_run_stats.inrange,
            "strand": minion_run_stats.strand,
            "pore": minion_run_stats.pore,
            "no_pore": minion_run_stats.no_pore,
            "zero": minion_run_stats.zero,
            "occupancy": minion_run_stats.occupancy(),
            "minKNOW_read_count": minion_run_stats.minKNOW_read_count,
            "minKNOW_histogram_values": minion_run_stats.minKNOW_histogram_values,
            "minKNOW_histogram_bin_width": minion_run_stats.minKNOW_histogram_bin_width,
        }

        result.append(element)

    myresult = dict()

    myresult["data"] = result

    if minion_run_status:
        myresult["minKNOW_colours_string"] = minion_run_status.minKNOW_colours_string

    else:
        myresult["minKNOW_colours_string"] = None

    return JsonResponse(myresult, safe=False)


@api_view(["GET"])
def flowcell_tasks_detail_all(request, pk):
    flowcell = Flowcell.objects.get(pk=pk)

    run_list = Run.objects.filter(flowcell=flowcell)

    queryset = JobType.objects.filter(private=False)

    result = []

    for jobtype in queryset:

        obj = {}
        obj.update(
            {
                "name": jobtype.name,
                "description": jobtype.description,
                "long_description": jobtype.long_description,
                "read_count": jobtype.readcount,
                "reference": jobtype.reference,
                "transcriptome": jobtype.transcriptome,
            }
        )

        jobmasterlist = JobMaster.objects.filter(
            Q(run__in=run_list) | Q(flowcell=flowcell)
        ).filter(job_type=jobtype)

        if len(jobmasterlist) > 0:
            obj2 = {}
            if jobmasterlist[0].reference:
                reference_name = jobmasterlist[0].reference.reference_name
            else:
                reference_name = ""
            obj2.update(
                {
                    "reference": reference_name,
                    "last_read": jobmasterlist[0].last_read,
                    "read_count": jobmasterlist[0].read_count,
                    "temp_file": jobmasterlist[0].tempfile_name,
                    "complete": jobmasterlist[0].complete,
                    "running": jobmasterlist[0].running,
                }
            )

            if len(obj2) > 0:
                obj.update({"job_details": obj2})

        result.append(obj)

    return HttpResponse(json.dumps(result), content_type="application/json")


@api_view(["GET"])
def flowcell_tabs_list(request, pk):
    """
    Get the flowcell tabs for a given flowcell.
    Parameters
    ----------
    request: rest_framework.request.Request
    pk: int
        Primary key of Flowcell record.

    Returns
    -------
    list
        List of tabs with data for flowcell.

    """

    tabs = ["summary-data", "tasks", "sharing", "notifications"]

    #
    # Check for basecalled data
    #
    flowcell_summary_barcode_list = FlowcellSummaryBarcode.objects.filter(
        flowcell__id=pk
    )

    if flowcell_summary_barcode_list.count() > 0:
        tabs.append("basecalled-data")
        tabs.append("reads")

    minion_run_stats_list = MinionRunStats.objects.filter(run_id__flowcell__id=pk)

    if minion_run_stats_list.count() > 0:
        tabs.append("live-event-data")

    paf_rough_cov_list = PafRoughCov.objects.filter(flowcell_id=pk)

    if paf_rough_cov_list.count() > 0:
        tabs.append("sequence-mapping")

    chromosome = ExpectedBenefitChromosomes.objects.filter(task__flowcell__id=pk).last()

    if chromosome is not None:
        tabs.append("advanced-sequence-mapping")

    centrifuge_output_list = CentrifugeOutput.objects.filter(task__flowcell_id=pk)

    if centrifuge_output_list.count() > 0:
        tabs.append("metagenomics")

    # Check for assembly data

    # TODO add an actual check here
    if JobMaster.objects.filter(
            flowcell__id=pk, job_type__name="Track Artic Coverage"
    ).last():
        if (
                JobMaster.objects.filter(
                    flowcell__id=pk, job_type__name="Track Artic Coverage"
                )
                        .last()
                        .read_count
        ):
            tabs.append("artic")

    return Response(tabs)


@api_view(["GET", "POST"])
def grouprun_list(request):
    """
    Not implemented
    Parameters
    ----------
    request

    Returns
    -------

    """

    if request.method == "GET":

        queryset = GroupRun.objects.distinct().filter(owner=request.user)

        serializer = GroupRunSerializer(
            queryset, many=True, context={"request": request}
        )

        return Response(serializer.data)

    elif request.method == "POST":

        serializer = GroupRunSerializer(data=request.data, context={"request": request})

        if serializer.is_valid():
            serializer.save(owner=request.user)

            return Response(serializer.data, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(["GET"])
def grouprun_detail(request, pk):
    """
    nOT IMPLEMENTED
    Parameters
    ----------
    request
    pk

    Returns
    -------

    """

    search_criteria = request.GET.get("search_criteria", "id")

    if search_criteria == "id":

        grouprun_list = GroupRun.objects.filter(owner=request.user).filter(id=pk)

    elif search_criteria == "name":

        grouprun_list = GroupRun.objects.filter(owner=request.user).filter(name=pk)

    else:

        grouprun_list = GroupRun.objects.none()

    if len(grouprun_list) != 1:
        return Response(status=status.HTTP_404_NOT_FOUND)

    grouprun = grouprun_list[0]

    serializer = GroupRunSerializer(grouprun, context={"request": request})

    return Response(serializer.data)


@api_view(["GET", "POST"])
def grouprun_membership_list(request):
    """
    Not implemented
    Parameters
    ----------
    request

    Returns
    -------

    """

    if request.method == "GET":

        grouprun_id = request.GET.get("grouprun_id")

        queryset = GroupRun.objects.filter(id=grouprun_id).filter(owner=request.user)

        if len(queryset) == 1:
            serializer = RunSerializer(
                queryset.first().runs, many=True, context={"request": request}
            )

            return Response(serializer.data)

    elif request.method == "POST":

        grouprun_id = request.data["grouprun_id"]
        run_id = request.data["run_id"]

        grouprun = GroupRun.objects.get(pk=grouprun_id)
        run = Run.objects.get(pk=run_id)

        run.groupruns.add(grouprun)
        run.save()

        serializer = RunSerializer(run, context={"request": request})

        return Response(serializer.data, status=status.HTTP_201_CREATED)


@api_view(["GET", "POST"])
def read_list_new(request):
    """
    API endpoint for either getting a list of reads to display in the read data table,
    or Posting new reads from the client.
    :param request: Django rest framework request object
    :type request: rest_framework.request.Request
    :return:
    """

    search_criteria = request.GET.get("search_criteria", "name")

    search_value = request.GET.get("search_value", "name")

    offset = int(request.GET.get("offset", "0"))

    limit = int(request.GET.get("limit", "10"))

    if request.method == "GET":

        if search_criteria == "run":

            qs = FastqRead.objects.filter(run__id=search_value)[offset: offset + limit]

        elif search_criteria == "fastqfile":

            qs = FastqRead.objects.filter(fastqfile_id=search_value)[
                 offset: offset + limit
                 ]

        else:

            qs = FastqRead.objects.none()

        serializer = FastqReadGetSerializer(qs, many=True, context={"request": request})

        return Response(serializer.data)

    elif request.method == "POST":

        reads = request.data
        if not reads:
            return Response("No reads in post request", status=status.HTTP_200_OK)
        save_reads_bulk(reads)
        return Response(
            "The task has started. It's in gods hands now...",
            status=status.HTTP_201_CREATED,
        )


@api_view(["GET"])
def readextra_list(request):
    search_criteria = request.GET.get("search_criteria", "read_id")

    search_value = request.GET.get("search_value", "name")

    if search_criteria == "read_id":

        qs = FastqRead.objects.filter(read_id=search_value)

        if len(qs) > 0:

            fastq_record = qs[0]

            fastq_string = ">{} runid={} ch={} read={} start_time={}\n{}\n".format(
                fastq_record.read_id,
                fastq_record.run.runid,
                fastq_record.channel,
                fastq_record.read,
                fastq_record.start_time.strftime("%Y-%m-%dT%H:%M:%SZ"),
                fastq_record.sequence,
            )

            read = {"sequence": fastq_string}

        else:
            read = {}

        return JsonResponse(read)


@api_view(["POST", "GET"])
def barcode_list_new(request):
    """

    Parameters
    ----------
    request

    Returns
    -------

    """

    if request.method == "GET":

        barcode_list = Barcode.objects.all()

        serializer = BarcodeSerializer(
            barcode_list, many=True, context={"request": request}
        )

        return Response(serializer.data)

    elif request.method == "POST":

        serializer = BarcodeSerializer(data=request.data, context={"request": request})

        if serializer.is_valid():
            serializer.save()

            return Response(serializer.data, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(["GET"])
def version(request):
    resp = {"server": "1.0", "clients": ["1.0"], "minknow": ["3.6.5"]}

    return HttpResponse(json.dumps(resp), content_type="application/json")


@api_view(["GET", "POST"])
def flowcell_sharing(request, pk):
    """
    list and add user permissions on a flowcell
    """

    flowcell = Flowcell.objects.get(pk=pk)

    if request.method == "POST":

        try:
            user = User.objects.get(username=request.data["username"])

        except User.DoesNotExist:

            return Response(
                {"message": "User does not exist."}, status=status.HTTP_404_NOT_FOUND
            )

        if flowcell.owner.id != request.user.id:
            return Response(
                {"message": "You are not authorised to share this flowcell."},
                status=status.HTTP_404_NOT_FOUND,
            )

        permission = request.data["permission"].lower()

        if user.has_perm(permission, flowcell):

            return Response(
                {
                    "message": "User already has permission. If you want to change the permission"
                               " level, delete the current permission first."
                },
                status=status.HTTP_404_NOT_FOUND,
            )

        else:

            assign_perm(permission, user, flowcell)

            return Response({}, status=status.HTTP_201_CREATED)

    else:

        users = User.objects.all()

        permission_list = []

        for user in users:

            if user != flowcell.owner:

                if user.has_perm("view_data", flowcell):
                    permission_list.append(
                        {
                            "username": user.username,
                            "permission": "VIEW DATA",
                            "permission_code": "view_data",
                            "user": user.id,
                            "flowcell": flowcell.id,
                        }
                    )

                if user.has_perm("run_analysis", flowcell):
                    permission_list.append(
                        {
                            "username": user.username,
                            "permission": "RUN ANALYSIS",
                            "permission_code": "run_analysis",
                            "user": user.id,
                            "flowcell": flowcell.id,
                        }
                    )

        return Response(permission_list)


@api_view(["POST"])
def flowcell_sharing_delete(request, pk):
    """
    delete a user from a flowcell list of sharing
    """

    flowcell = Flowcell.objects.get(pk=pk)

    data = request.data

    try:
        user = User.objects.get(pk=data["user_id"])

    except User.DoesNotExist:

        return Response(
            {"message": "User does not exist."}, status=status.HTTP_404_NOT_FOUND
        )

    permission = data["permission"]

    if request.user == user and user.has_perm(permission, flowcell):

        remove_perm(permission, user, flowcell)
        return Response({"message": "Permission deleted"}, status=status.HTTP_200_OK)

    elif request.user == flowcell.owner:

        remove_perm(permission, user, flowcell)
        return Response({"message": "Permission deleted"}, status=status.HTTP_200_OK)

    else:

        return Response(
            {"message": "You do not have the permission to execute this action."},
            status=status.HTTP_400_BAD_REQUEST,
        )


@api_view(["GET", "POST"])
def job_master_list(request):
    """
    API endpoint for dealing with fetching or creating new JobMasters
    :param request: The django rest framework request body - if a GET request should contain
    search_criteria and search_value params
    :return: A dictionary containing all the matching JobMasters on this flowcell
    """
    # TODO this is hacky as well
    if request.method == "GET":
        search_criteria = request.GET.get("search_criteria", "flowcell")
        if search_criteria == "flowcell":
            flowcell_id = request.GET.get("search_value", -1)
            # If there was no flowcell ID provided
            if flowcell_id == -1:
                return Response(
                    "No flowcell ID provided", status=status.HTTP_400_BAD_REQUEST
                )
            # Get all the JobMasters for this flowcell, excluding job type other - which are for background sub tasks
            tasks_list = JobMaster.objects.filter(
                flowcell__id=int(flowcell_id)
            ).exclude(job_type__name__in=["Other", "Run Artic"])
            # Serialise the data to a python object
            serializer = JobMasterSerializer(tasks_list, many=True)
            # Check the icon required on the pause control part of the table
            for data in serializer.data:
                if data["paused"]:
                    data["icon"] = "play"
                    data["iconText"] = "Play"
                else:
                    data["icon"] = "pause"
                    data["iconText"] = "Pause"
            result = {"data": serializer.data}
            return Response(result, status=200)
        else:
            return Response(
                "Bad search criteria provided", status=status.HTTP_400_BAD_REQUEST
            )
    # Its a post request to create a new task
    else:
        try:
            # Get the flowcell by name
            flowcell = Flowcell.objects.get(
                Q(name=request.data["flowcell"])
                | Q(pk=int(request.data["flowcell"]))
            )
            request.data["flowcell"] = flowcell.id
        # If that doesn't work
        except Flowcell.DoesNotExist as e:
            print("Exception: {}".format(e))
            return Response(
                "Flowcell not found. Please contact server admin.", status=status.HTTP_404_NOT_FOUND
            )
        if (
                not request.user.has_perm("run_analysis", flowcell)
                and not request.user == flowcell.owner
        ):
            return Response(
                {
                    "message": "Permission denied - you do not have the RUN ANALYSIS permission."
                },
                status=status.HTTP_401_UNAUTHORIZED,
            )
        # Check to see if we have a string as the flowcell, not an Int for a PK based lookup
        # This is for starting a job from the client

        # If the job isn't EB or minimap2 or artic
        if int(request.data["job_type"]) not in [4, 15, 16]:
            request.data["reference"] = None
            # Hard set the covid reference
        # if int(request.data["job_type"]) == 16:
        #     print("we are setting covid reference")
        #     reference = ReferenceInfo.objects.get(name="covid_19")
        #     request.data["reference"] = reference.id

        # For client side job starting, check if we have a reference and it's a string
        if (
            "reference" in request.data.keys()
            and type(request.data["reference"]) is str
        ):
            try:
                # If a reference hasn't been selected.
                if request.data["reference"] == "":
                    return Response(
                        "Reference does not exist - Please select a reference",
                        status=400,
                    )
                reference = ReferenceInfo.objects.get(
                    Q(name=request.data["reference"])
                    | Q(pk=int(request.data["reference"]))
                )
                request.data["reference"] = reference.id
            except ReferenceInfo.DoesNotExist as e:
                print("Exception: {}".format(e))
                return Response(
                    "Reference not found Please contact server admin", status=500
                )
        serializer = JobMasterInsertSerializer(data=request.data)
        # If the serialiser is valid
        if serializer.is_valid():
            task = serializer.save()
            response_data = {"message": "Task created successfully!", "pk": task.id}
            return JsonResponse(response_data)

        else:
            print(serializer.errors)
            return JsonResponse({"error_messages": serializer.errors}, status=500)


@api_view(["GET"])
def task_types_list(request):
    """
    Get the list of task types to display for Job creation options
    :param request: The request body object
    :return: A list of dictionaries containing the task id, name and description.
    """

    # If it's a request from the client
    if request.GET.get("cli", False):
        # These are the available tasks
        tasks = ["Metagenomics", "Minimap2", "Track Artic Coverage"]
        # Get the tasks
        result = JobType.objects.filter(name__in=tasks).values("id", "name", "description")

    else:
        # Otherwise it's for the site, so make all that are public available
        result = JobType.objects.filter(private=False).values("id", "name", "description")

    return Response({"data": result}, status=status.HTTP_200_OK)


@api_view(["POST"])
def task_control(request):
    """
    Endpoint for clearing a task and resetting it to start form the beginning
    :param request: Request object, contains the Id of the task keyed to flowcellJobId
    :return: A status reflecting the state of the executed code
    """
    # Get the task object from django ORM

    job_master = JobMaster.objects.get(pk=request.data["flowcellJobId"])
    if not (
            request.user == job_master.flowcell.owner
            or request.user.has_perm("run_analysis", job_master.flowcell)
    ):
        return Response(
            "You do not have permission to perform this action.", status=403
        )

    action_type = request.data["actionType"]
    lookup_action_type = {1: "Reset", 2: "Pause", 3: "Delete"}
    action = lookup_action_type[action_type]
    unrecognised_action_message = (
        "Apologies, but this action type was not recognised."
        "It may have not been implemented yet."
    )

    if job_master.job_type.name == "Metagenomics":
        if action == "Reset":
            delete_metagenomics_task.delay(job_master.id, True)
            return_message = (
                f"Successfully started to reset this metagenomics task, id: {job_master.id}."
                f" Clearing previous data may take a while, please be patient!"
            )
        elif action == "Pause":
            job_master, return_message = pause_job(job_master)
        elif action == "Delete":
            delete_metagenomics_task.delay(job_master.id)
            return_message = (
                f"Successfully started deletion of this metagenomics task, id: {job_master.id}."
                f" Clearing previous data may take a while, please be patient!"
            )
        else:
            return Response(unrecognised_action_message, status=500)

    elif job_master.job_type.name == "Minimap2":
        if action == "Reset":
            delete_alignment_task.delay(job_master.id, True)
            return_message = (
                f"Successfully began reset of minimap2 task, id: {job_master.id}"
                f" Clearing previous data may take a while, please be patient!"
            )
        elif action == "Pause":
            job_master, return_message = pause_job(job_master)
        elif action == "Delete":
            delete_alignment_task.delay(job_master.id)
            return_message = (
                f"Successfully began deletion of minimap2 task, id: {job_master.id}"
                f" Clearing previous data may take a while, please be patient!"
            )
        else:
            return Response(unrecognised_action_message, status=500)
    # for artic tasks
    elif job_master.job_type.id == 16:
        if action == "Reset":
            clear_artic_data(job_master)
            clear_artic_command_job_masters(job_master.flowcell_id)
            ArticBarcodeMetadata.objects.filter(job_master=job_master).delete()
            job_master.read_count = 0
            job_master.last_read = 0
            job_master.running = False
            job_master.complete = False
            job_master.save()
            return_message = "Successfully reset artic Task."
        elif action == "Pause":
            job_master, return_message = pause_job(job_master)
        elif action == "Delete":
            clear_artic_command_job_masters(job_master.flowcell_id)
            clear_artic_data(job_master)
            job_master.delete()
            return_message = "Successfully deleted artic Task."

    else:
        return Response("Not implemented, please deal with as an admin.", status=status.HTTP_501_NOT_IMPLEMENTED)

    return Response(return_message, status=200)
