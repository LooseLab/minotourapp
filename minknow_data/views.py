import datetime
import json
from collections import defaultdict

import pytz
from django.http import HttpResponse
# Create your views here.
from django.utils import timezone
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from minknow_data.models import Minion, MinionEvent, MinionRunStats, MinionControl, MinionScripts, MinionInfo, \
    MinionMessage, MinionRunInfo, Run, Flowcell, MinionEventType
from minknow_data.serializers import MinionRunStatsSerializer, MinionSerializer, MinIONControlSerializer, \
    MinIONEventSerializer, MinIONScriptsSerializer, MinionInfoSerializer, MinionMessageSerializer, \
    MinionRunInfoSerializer, MinIONEventTypeSerializer
from reads.utils import create_histogram_series


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
                    "ACQUISITION_COMPLETED",
                    "ACQUISITION_STARTING",
                    "No Run",
                    "ACQUISITION_RUNNING",
                    "ACQUISITION_PROCESSING",
                    "ACQUISITION_FINISHING",
                    "Nothing running"
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
            data_to_add = data.get(active_minion["name"], blank_stats)
            # Remove the minion run stats id so we don't overwrite the minion ID
            data_to_add.pop("id")
            active_minion.update(data.get(active_minion["name"], blank_stats))
            active_minion.update(extra_data[active_minion["name"]])

    return Response(serializer.data)


@api_view(["GET", "POST"],)
def minion_control_list(request, pk):
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


@api_view(["GET", "POST"],)
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
    ["GET"]
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


@api_view(["GET", "PUT", "POST", "DELETE"])
def minion_info_list(request, pk):
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


@api_view(["GET", "POST"])
def minION_detail(request, pk):
    # """
    # TODO describe function
    # """

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
        serializer = MinionSerializer(data=request.data, context={"request": request})
        if serializer.is_valid():
            serializer.save(owner=request.user)
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


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


@api_view(["GET"],)
def recentminion_messages_list(request, pk):
    queryset = MinionMessage.objects.filter(minION=pk).filter(
        minKNOW_message_timestamp__gte=timezone.now() - datetime.timedelta(hours=24)
    )

    serializer = MinionMessageSerializer(
        queryset, many=True, context={"request": request}
    )

    return Response(serializer.data)


@api_view(["GET"])
def minknow_message(request, pk):
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
    return_list = []
    for message in messages:
        message_text = message.message if not message.full_text else message.full_text
        return_list.append(
            {"messageText": message_text, "timestamp": message.timestamp}
        )
    return Response({"data": return_list})


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
def minknow_control_run_stats(request):
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
def minion_run_status_list(request, pk):
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
    if request.method == "POST":
        serializer = MinionRunInfoSerializer(
            data=request.data, context={"request": request}
        )
        if serializer.is_valid():
            serializer.save()
            # update the run start time
            run = Run.objects.get(pk=pk)
            try:
                run_start_date = datetime.datetime.strptime(serializer.data["minKNOW_start_time"],
                                                    "%Y-%m-%dT%H:%M:%S%z")
            except ValueError as e:
                run_start_date = datetime.datetime.strptime(
                    serializer.data["minKNOW_start_time"], "%Y-%m-%dT%H:%M:%S"
                ).replace(tzinfo=pytz.UTC)
            if run_start_date < run.start_time:
                run.start_time = run_start_date
                run.save()
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


@api_view(["GET"])
def test_connection(request):
    """
    Test a connection from the minFQ client
    Parameters
    ----------
    request: rest_framework.request.Request
        Ajax request object
    Returns
    -------
    rest_framework.response.Response
        A status code of 200.
    """
    return Response("Hello from the other siddde. Here we can use F strings. Wow", status=status.HTTP_200_OK)


@api_view(["GET"])
def version(request):
    """
    Return the version of the server being run, and the supported version of minknow and minFQ
    Parameters
    ----------
    request: rest_framework.request.Request
        The request
    Returns
    -------
    json
        Dict of server version, clients supported and minknow version supported
    """
    resp = {"server": "1.0", "clients": ["1.0", "1.1", "1.2", "1.3"], "minknow": ["3.6.5", "4.0.3", "4.0.4", "4.0.5"]}
    return HttpResponse(json.dumps(resp), content_type="application/json")


@api_view(["GET"])
def flowcell_minknow_stats_list(request, pk, check_id):
    """
    Fetch pre parsed minknow metrics about this flowcell
    Parameters
    ----------
    request: rest_framework.request.Request
        The ajax request body
    pk: int
        Primary key of the flowcell
    check_id: int
        Last read id to read from.

    Returns
    -------
    dict
        Dictionary of prepared minKnow stats to be displayed on the live data page
    """
    # Take the first entry to be the run status
    minion_run_meta_information = (
        MinionRunInfo.objects.filter(run_id__flowcell_id=pk)
        .last()
    )
    flowcell = Flowcell.objects.get(pk=pk)
    runs = Run.objects.filter(flowcell=flowcell)
    run_info = [(name, (run_start.timestamp() * 1000)) for name, run_start in runs.values_list("name", "start_time") if run_start is not None]
    if not MinionRunStats.objects.filter(
        run_id__in=flowcell.runs.all().exclude(name="mux_scan"), id__gt=check_id
    ).count():
        return Response({"last_minion_run_stats_id": check_id}, status=status.HTTP_206_PARTIAL_CONTENT)
    minion_run_stats_gen = (
        x
        for x in list(
            MinionRunStats.objects.filter(
                run_id__in=flowcell.runs.all().exclude(name="mux_scan"), id__gt=check_id
            )
        )
    )
    result = defaultdict(list)
    result_pore = defaultdict(dict)
    result_strand = defaultdict(dict)
    result_temperature = defaultdict(dict)
    possible_pore_states = {
        "above",
        "adapter",
        "below",
        "good_single",
        "strand",
        "inrange",
        "multiple",
        "pending_mux_change",
        "saturated",
        "unavailable",
        "unblocking",
        "unclassified",
        "unknown",
        "zero",
        "pore",
        "no_pore",
    }
    strand_counts = ["strand", "pore", "adapter", "good_single"]
    result_temperature["temperature_history"]["asic_temp"] = {
        "name": "Asic temp",
        "data": [],
    }
    result_temperature["temperature_history"]["heat_sink_temp"] = {
        "name": "HeatSink temp",
        "data": [],
    }
    minknow_colours = {}
    # Create a lookup dictionary of labels and styles for each pore state
    for x in json.loads(minion_run_meta_information.minKNOW_colours_string)["groups"]:
        if "states" in x:
            for j in x["states"]:
                minknow_colours[j["name"]] = j.get(
                    "style", {"label": j["name"], "colour": "000000"}
                )
        else:
            minknow_colours[x] = {"label": j["name"], "colour": "000000"}
    # Manually add in some of the other that aren't recorded anymore, should be fixed later
    minknow_colours["good_single"] = {"label": "Good single", "colour": "32CD32"}
    # Create the entry in the results dicionary for each pore state
    # TODO check the pore states match
    last_time = MinionRunStats.objects.filter(run__flowcell_id=pk).last().sample_time
    first_time = MinionRunStats.objects.filter(run__flowcell_id=pk).first().sample_time
    div, rem = divmod((last_time - first_time).total_seconds(), 24*3600)
    counter = div
    for x in possible_pore_states:
        _ = {
            "name": minknow_colours.get(x, {"label": x})["label"],
            "color": f"#{minknow_colours.get(x, {'colour': '000000'})['colour']}",
            "data": [],
        }
        result_pore[x] = _
        if x in strand_counts:
            result_strand[x] = _

    histogram_series_suffix = (
        "Bases" if minion_run_meta_information.local_basecalling else "ev"
    )
    # Loop our results
    for mrs in minion_run_stats_gen:
        if not div == 0 and counter % div == 0:
            counter += 1
            continue
        else:
            counter += 1
        sample_time_microseconds = int(
            mrs.sample_time.replace(microsecond=0).timestamp() * 1000
        )
        result["yield_history"].append([sample_time_microseconds, mrs.event_yield])
        for x in possible_pore_states:
            a = [sample_time_microseconds, mrs.__dict__.get(x, 0)]
            result_pore[x]["data"].append(a)
            # strand calculations
            # if x in strand_counts:
            #     result_strand[x]["data"].append(a)
        if mrs.minKNOW_histogram_values:
            histogram_values = (
                mrs.minKNOW_histogram_values.replace("'", "")
                .replace('"', "")[1:-1]
                .split(", ")
            )
            histogram_bin_width = int(mrs.minKNOW_histogram_bin_width)
            n50 = int(mrs.n50_data)
            result["histogram_history"].append(
                {
                    "histogram_series": [
                        create_histogram_series(
                            int(value),
                            n50,
                            histogram_bin_width,
                            ind,
                            histogram_series_suffix,
                        )
                        for ind, value in enumerate(histogram_values)
                    ],
                    "histogram_bin_with": histogram_bin_width,
                    "n50_data": n50,
                    "sample_time": mrs.sample_time,
                    "categories": [
                        f"{x*histogram_bin_width} - {(x+1)*histogram_bin_width} ev"
                        for x in range(len(histogram_values))
                    ],
                }
            )
        result_temperature["temperature_history"]["asic_temp"]["data"].append(
            [sample_time_microseconds, round(mrs.asic_temp, 2)]
        )
        result_temperature["temperature_history"]["heat_sink_temp"]["data"].append(
            [sample_time_microseconds, round(mrs.heat_sink_temp)]
        )
        result["voltage_history"].append([sample_time_microseconds, mrs.voltage_value])
        result["occupancy_history"].append([sample_time_microseconds, mrs.occupancy()])
        last_minion_run_stats_id = mrs.id
    # TOdo here is where we need a check for events or basecalled data, maybe add to model
    yield_name = "Bases" if minion_run_meta_information.local_basecalling else "Events"
    result["yield_history"] = [{"name": yield_name, "data": result["yield_history"]}]
    result["voltage_history"] = [{"name": "Voltage", "data": result["voltage_history"]}]
    result["occupancy_history"] = [
        {"name": "% Occupancy", "data": result["occupancy_history"]}
    ]
    result["pore_history"] = result_pore.values()
    result["in_strand_history"] = result_strand.values()
    result["temperature_history"] = result_temperature["temperature_history"].values()
    result["last_minion_run_stats_id"] = last_minion_run_stats_id
    result["run_info"] = run_info
    return Response(result, status=200)


@api_view(["GET"])
def events_type_list(request):
    """
    :purpose: Returns a list of all possible events to the client
    :used_by: minotour client

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


