import datetime
import json
from datetime import timedelta

import dateutil.parser
from celery import task
from celery.utils.log import get_task_logger
from dateutil import parser
from django.core.paginator import EmptyPage, PageNotAnInteger, Paginator
from django.core.serializers.json import DjangoJSONEncoder
from django.db.models import Max
from django.db.models import Q
from django.http import HttpResponse, JsonResponse
from django.shortcuts import render
from django.utils import timezone
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.views import APIView

from alignment.models import PafRoughCov
from centrifuge.models import CentrifugeOutput
from jobs.models import JobMaster, JobType
from minotourapp import settings
from reads.models import (Barcode, FastqFile, FastqRead, FastqReadType,
                          MinIONControl, MinIONEvent,
                          MinIONEventType, MinionMessage, MinIONRunStats,
                          MinIONRunStatus, MinIONScripts, MinIONStatus, Run, GroupRun, FlowcellStatisticBarcode,
                          FlowcellSummaryBarcode, Flowcell, MinION, FlowcellTab, RunSummary)
from reads.models import FlowcellChannelSummary
from reads.models import FlowcellHistogramSummary
from reads.serializers import (BarcodeSerializer,
                               ChannelSummarySerializer, FastqFileSerializer, FastqReadSerializer,
                               FastqReadTypeSerializer, FlowcellSerializer, MinIONControlSerializer,
                               MinIONEventSerializer,
                               MinIONEventTypeSerializer,
                               MinionMessageSerializer,
                               MinIONRunStatsSerializer,
                               MinIONRunStatusSerializer,
                               MinIONScriptsSerializer, MinIONSerializer,
                               MinIONStatusSerializer,
                               RunSerializer,
                               RunStatisticBarcodeSerializer,
                               RunSummaryBarcodeSerializer, ChannelSummary, RunStatisticBarcode, RunSummaryBarcode,
                               GroupRunSerializer, FlowcellSummaryBarcodeSerializer, FastqReadGetSerializer,
                               FlowcellTabSerializer)
from reads.utils import get_coords

logger = get_task_logger(__name__)

@api_view(['GET'])
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
    serializer = FastqReadTypeSerializer(queryset, many=True, context={'request': request})
    return Response(serializer.data)


@api_view(['GET'])
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

    serializer = FastqReadTypeSerializer(run, context={'request': request})
    return Response(serializer.data)


@api_view(['GET'])
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
    queryset = MinIONEventType.objects.all()
    serializer = MinIONEventTypeSerializer(queryset, many=True, context={'request': request})
    return Response(serializer.data)


@api_view(['GET'])
def events_type_detail(request, pk): # TODO consider removing
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
        event_ = MinIONEventType.objects.get(pk=pk)
    except MinIONEventType.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    serializer = MinIONEventTypeSerializer(event_, context={'request': request})
    return Response(serializer.data)



@api_view(['GET'])
def fastq_detail(request,pk):
    """

    :param request:
    :param pk: fastqid
    :return:
    """
    queryset = FastqFile.objects.filter(runid=pk)
    serializer = FastqFileSerializer(queryset, many=True, context={'request': request})
    return Response(serializer.data)


@api_view(['GET', 'POST'])
def fastq_file(request,pk):
    """
    :purpose: returns an md5 checksum for a file as seen by minotour for a specific run id
    :used_by: minotour client gets data and checks observed files against it.
    :author: Matt Loose

    :param request: (standard django request) without querystring parameter
    :param pk: pk is the runid

    :return: (str) json format
    """
    '''if request.method == "GET":
        try:
            queryset = FastqFile.objects.filter(runid=pk)
            serializer = FastqFileSerializer(queryset, many=True, context={'request': request})
            return Response(serializer.data)
        except FastqFile.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)
    '''

    if request.method == 'GET':

        queryset = FastqFile.objects.filter(runid=pk).filter(owner=request.user)

        serializer = FastqFileSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)

    elif request.method == 'POST':

        run = Run.objects.filter(runid=request.data["runid"]).filter(owner=request.user).first()

        obj,created = FastqFile.objects.get_or_create(
            runid = request.data["runid"],
            name = request.data["name"],
            owner = request.user,
        )
        obj.run = run
        obj.md5 = request.data["md5"]


        obj.save()

        serializer = FastqFileSerializer(obj, context={'request': request})

        return Response(serializer.data, status=status.HTTP_201_CREATED)

    return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET', 'POST'])
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

    if request.method == 'GET':

        queryset = Run.objects.filter(owner=request.user).filter(to_delete=False)
        serializer = RunSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)

    elif request.method == 'POST':

        serializer = RunSerializer(
            data=request.data,
            context={'request': request}
        )

        if serializer.is_valid():

            serializer.save(owner=request.user)
            return Response(serializer.data, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET', 'PUT'])
def run_detail(request, pk):
    """
    Retrieve, update or delete a run instance.
    """

    search_criteria = request.GET.get('search_criteria', 'id')

    if search_criteria == 'runid':

        run_list = Run.objects.filter(owner=request.user).filter(runid=pk)

    elif search_criteria == 'id':

        run_list = Run.objects.filter(owner=request.user).filter(id=pk)

    else:

        run_list = Run.objects.none()

    if len(run_list) < 1:

        return Response(status=status.HTTP_404_NOT_FOUND)

    run = run_list[0]  # TODO what if by any means we have more than one run if the same runid?

    if request.method == 'GET':
        serializer = RunSerializer(run, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'PUT':
        serializer = RunSerializer(run, data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET', 'POST'])
def minion_list(request):
    """
    List of all minIONs by user, or create a new minION.
    """
    if request.method == 'GET':
        queryset = MinION.objects.filter(owner=request.user)
        serializer = MinIONSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'POST':
        serializer = MinIONSerializer(data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save(owner=request.user)
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET'])
def activeminion_list(request):
    """
    List of all minIONs by user, or create a new minION.
    """

    active_minion_list = []

    for minion in MinION.objects.filter(owner=request.user):

        minion_event_list = MinIONEvent.objects.filter(minION=minion).order_by('datetime')

        if minion_event_list.count() > 0:

            last_minion_event = minion_event_list.last()

            if last_minion_event.event.name != 'unplugged':

                active_minion_list.append(minion)

    serializer = MinIONSerializer(active_minion_list, many=True, context={'request': request})

    return Response(serializer.data)


@api_view(['GET', 'POST'])
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

    if request.method == 'GET':
        queryset = MinionMessage.objects.filter(minION=pk)
        serializer = MinionMessageSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'POST':
        serializer = MinionMessageSerializer(data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET'],)
def sinceminion_messages_list(request, pk, starttime, endtime):

    correctedstart = parser.parse(starttime) - timedelta(minutes=180)

    correctedend = parser.parse(endtime) + timedelta(minutes=180)

    queryset = MinionMessage.objects.filter(minION=pk)\
        .filter(minKNOW_message_timestamp__gte=correctedstart.isoformat().replace('+00:00', 'Z'))\
        .filter(minKNOW_message_timestamp__lte=correctedend.isoformat().replace('+00:00', 'Z'))

    serializer = MinionMessageSerializer(queryset, many=True, context={'request': request})

    return Response(serializer.data)


@api_view(['GET'],)
def recentminion_messages_list(request, pk):

    queryset = MinionMessage.objects.filter(minION=pk)\
        .filter(minKNOW_message_timestamp__gte=timezone.now() - timedelta(hours=24))

    serializer = MinionMessageSerializer(queryset, many=True, context={'request': request})

    return Response(serializer.data)


@api_view(['GET'])
def minknow_message_list_by_flowcell(request, pk):
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
    form_start_time = request.GET.get('start_time', None)
    form_end_time = request.GET.get('end_time', None)

    start_time = None
    end_time = None

    if form_start_time:

        try:

            start_time = dateutil.parser.parse(form_start_time)

        except ValueError:

            print('Error when parsing start_time')

    if form_end_time:

        try:

            end_time = dateutil.parser.parse(form_end_time)

        except ValueError:

            print('Error when parsing end_time')

    flowcell = Flowcell.objects.get(pk=pk)

    minion_list = []

    for run in flowcell.runs.all():

        if run.minion:

            minion_list.append(run.minion)

    if start_time and end_time:

        messages = MinionMessage.objects\
            .filter(minion__in=minion_list)\
            .filter(timestamp__gt=start_time)\
            .filter(timestamp__lt=end_time) \
            .order_by('-timestamp')

    elif start_time:

        messages = MinionMessage.objects\
            .filter(minion__in=minion_list)\
            .filter(timestamp__gt=start_time) \
            .order_by('-timestamp')

    elif end_time:

        messages = MinionMessage.objects\
            .filter(minion__in=minion_list)\
            .filter(timestamp__lt=end_time) \
            .order_by('-timestamp')

    else:

        messages = MinionMessage.objects\
            .filter(minion__in=minion_list) \
            .order_by('-timestamp')[:10]

    return render(request, 'reads/minknow_messages.html', {'message_list': messages})


@api_view(['GET','POST'],)
def minION_control_list(request,pk):
    """
    TODO describe function
    """
    if request.method == 'GET':
        queryset = MinIONControl.objects.filter(minION=pk).filter(complete=False)
        serializer = MinIONControlSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'POST':
        serializer = MinIONControlSerializer(data=request.data, context={'request': request})

        if serializer.is_valid():
            serializer.save(owner=request.user)
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET','POST'],)
def minION_control_update(request, pk, checkid):
    """
    TODO describe function
    """
    try:
        event_ = MinIONControl.objects.get(id=checkid)

        if event_.owner != request.user:
            return Response(status=status.HTTP_401_UNAUTHORIZED)

    except MinIONControl.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        serializer = MinIONControlSerializer(event_, context={'request': request})
        return Response(serializer.data)

    if request.method == 'POST':

        event_ = MinIONControl.objects.get(id=checkid)
        event_.complete=True
        event_.save(update_fields=['complete'])
        return Response(status=status.HTTP_204_NO_CONTENT)


@api_view(['GET',])
def minION_currentrun_list(request, pk):
    """
    TODO describe function
    """
    try:
        minion = MinION.objects.get(pk=pk)
    except MinION.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        serializer = MinIONSerializer(minion, context={'request': request})
        return Response(serializer.data)


@api_view(['GET', 'POST'])
def minION_events_list(request, pk):
    """
    List of all events associated with a minION or add a new event.
    """
    if request.method == 'GET':
        queryset = MinIONEvent.objects.filter(minION=pk)
        serializer = MinIONEventSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'POST':
        serializer = MinIONEventSerializer(data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET', 'POST'])
def minION_scripts_list(request, pk):
    """
    List of all scripts associated with a minION.
    """
    if request.method == 'GET':
        queryset = MinIONScripts.objects.filter(minION=pk)
        serializer = MinIONScriptsSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'POST':
        serializer = MinIONScriptsSerializer(data=request.data, partial=True, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET', 'PUT', 'POST', 'DELETE'])
def minION_status_list(request,pk):
    """
    TODO describe function
    """
    if request.method == 'POST':
        serializer = MinIONStatusSerializer(data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    try:
        minIONstat = MinIONStatus.objects.get(minION=pk)

    except MinIONStatus.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        serializer = MinIONStatusSerializer(minIONstat, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'PUT':
        serializer = MinIONStatusSerializer(minIONstat, data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    elif request.method == 'DELETE':
        minIONstat.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


@api_view(['GET', 'PUT', 'POST', 'DELETE'])
def minION_run_stats_list(request,pk):
    """
    TODO describe function
    """
    if request.method == 'POST':
        serializer = MinIONRunStatsSerializer(data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    try:
        crazyminIONrunstats = MinIONRunStats.objects.filter(run_id=pk)

    except MinIONRunStats.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':

        serializer = MinIONRunStatsSerializer(crazyminIONrunstats, many=True , context={'request': request})

        return Response(serializer.data)

    elif request.method == 'PUT':
        serializer = MinIONRunStatsSerializer(crazyminIONrunstats, data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    elif request.method == 'DELETE':
        crazyminIONrunstats.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


@api_view(['GET'])
def minION_run_stats_latest(request,pk,checkid):
    """
    TODO describe function
    """
    try:
        crazyminIONrunstats = MinIONRunStats.objects.filter(run_id=pk, id__gt=checkid)[:1000]

    except MinIONRunStats.DoesNotExist:

        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':

        serializer = MinIONRunStatsSerializer(crazyminIONrunstats, many=True , context={'request': request})

        return Response(serializer.data)


@api_view(['GET', 'PUT', 'POST', 'DELETE'])
def minION_run_status_list(request,pk):
    """
    TODO describe function
    """
    if request.method == 'POST':

        serializer = MinIONRunStatusSerializer(data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    try:
        minIONrunstat = MinIONRunStatus.objects.filter(run_id=pk)

    except MinIONRunStatus.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        serializer = MinIONRunStatusSerializer(minIONrunstat, many=True, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'PUT':
        serializer = MinIONRunStatusSerializer(minIONrunstat, data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    elif request.method == 'DELETE':
        minIONrunstat.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


@api_view(['GET', 'PUT', 'DELETE'])
def minION_scripts_detail(request, pk, nk):
    """
    TODO describe function
    """
    try:
        script = MinIONScripts.objects.get(pk=nk)
    except MinIONScripts.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        serializer = MinIONScriptsSerializer(script, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'PUT':
        serializer = MinIONScriptsSerializer(script, data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    elif request.method == 'DELETE':
        script.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


@api_view(['GET', 'POST'])
def read_list(request, pk):

    if request.method == 'GET':

        querysets = FastqRead.objects.filter(run_id=pk)

        serializer = FastqReadGetSerializer(querysets, many=True, context={'request': request})

        return Response(serializer.data)

    elif request.method == 'POST':

        serializer = FastqReadSerializer(data=request.data, many=True, context={'request': request})

        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET', 'POST'])
def minION_detail(request, pk):
    # """
    # TODO describe function
    # """

    if request.method == 'GET':

        search_criteria = request.GET.get('search_criteria', 'id')

        if search_criteria == 'id':

            minion_list = MinION.objects.filter(owner=request.user).filter(id=pk)

        elif search_criteria == 'name':

            minion_list = MinION.objects.filter(owner=request.user).filter(name=pk)

        else:

            minion_list = MinION.objects.none()

        if len(minion_list) != 1:
            return Response(status=status.HTTP_404_NOT_FOUND)

        minion = minion_list[0]

        serializer = MinIONSerializer(minion, context={'request': request})

        return Response(serializer.data)

    elif request.method == 'POST':

        minion = MinION.objects.get(pk=pk)

        serializer = MinIONSerializer(minion, data=request.data, partial=True,context={'request': request})

        if serializer.is_valid():
            serializer.save()

        else:
            print ("PROBLEM")

        return Response(serializer.data)


@api_view(['GET'])
def readname_list(request, pk):
    """
    TODO describe function
    CUrrent API behaviour asssumes starting from page1.
    Request for page 0 returns last page.
    This could be a source of confusion and we should resolve.
    """
    if request.method == 'GET':
        queryset = FastqRead.objects.filter(fastqfile_id=pk).order_by('id')

        paginator = Paginator(queryset, settings.PAGINATION_PAGE_SIZE)

        page = request.GET.get('page')

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
        result3['number_pages'] = paginator.num_pages
        result3['data'] = list(result2)
        return HttpResponse(json.dumps(result3), content_type="application/json")


@api_view(['GET'])
def run_summary_barcode(request, pk):
    """
    Return a list with summaries for each read type of a given run.
    """
    queryset = RunSummaryBarcode.objects\
        .filter(run_id__owner=request.user)\
        .filter(run_id=pk)

    serializer = RunSummaryBarcodeSerializer(queryset, many=True, context={'request': request})

    return Response(serializer.data)


@api_view(['GET'])
def run_summary_barcode_by_minute(request, pk):
    """
    Return a list with summaries for a particular run grouped by minute.
    """
    queryset = RunStatisticBarcode.objects\
        .filter(run_id__owner=request.user)\
        .filter(run_id=pk)

    serializer = RunStatisticBarcodeSerializer(queryset, many=True, context={'request': request})

    return Response(serializer.data)


@api_view(['GET'])
def run_channel_summary(request, pk):
    """
    Return a list with channel info for a particular run.
    """
    queryset = ChannelSummary.objects\
        .filter(run_id__owner=request.user)\
        .filter(run_id=pk)

    serializer = ChannelSummarySerializer(queryset, many=True, context={'request': request})

    return Response(serializer.data)


@api_view(['GET', 'POST'])
def barcode_list(request, pk):
    """
    Return a list of barcodes of a particular run.
    """

    if request.method == 'GET':
        queryset = Barcode.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk)

        serializer = BarcodeSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)

    elif request.method == 'POST':

        barcodename=request.data["name"]

        minionrun = Run.objects.get(pk=pk)

        barcode, created2 = Barcode.objects.get_or_create(run=minionrun, name=barcodename)

        serializer = BarcodeSerializer(barcode, context={'request': request})

        return Response(serializer.data, status=status.HTTP_201_CREATED)


@api_view(['GET'])
def barcode_detail(request, pk):
    """
    Return the details of a particular barcode.
    """
    queryset = Barcode.objects.get(pk=pk)

    serializer = BarcodeSerializer(queryset, many=False, context={'request': request})

    return Response(serializer.data)


@api_view(['GET', 'POST'])
def flowcell_list(request):

    if request.method == 'GET':

        queryset = Flowcell.objects.filter(owner=request.user)

        flowcells = []

        for record in queryset:
            flowcell = {
                'id': record.id,
                'name': record.name,
                'size': record.size,
                'start_time': record.start_time,
                'number_reads': record.number_reads,
                'number_reads_processed': record.number_reads_processed,
                'number_runs': record.number_runs,
                'number_barcodes': record.number_barcodes,
                'total_read_length': record.total_read_length,
                'average_read_length': record.average_read_length,
                'is_active': record.is_active,
                'sample_name': record.sample_name,
                'has_fastq': record.has_fastq,
            }

            flowcells.append(flowcell)

        return JsonResponse({'data': flowcells})

    else:

        serializer = FlowcellSerializer(data=request.data, context={'request': request})

        if serializer.is_valid():

            serializer.save(owner=request.user)

            return Response(serializer.data, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET', 'POST'])
def flowcell_detail(request, pk):

    if request.method == 'GET':

        search_criteria = request.GET.get('search_criteria', 'id')

        if search_criteria == 'id':

            flowcell_list = Flowcell.objects.filter(owner=request.user).filter(id=pk)

        elif search_criteria == 'name':

            flowcell_list = Flowcell.objects.filter(owner=request.user).filter(name=pk)

        else:

            flowcell_list = Flowcell.objects.none()

        if len(flowcell_list) != 1:

            return Response({'data': {}})


        # TODO updated this, check with Roberto that this is cool

        flowcell = flowcell_list[0]

        serializer = FlowcellSerializer(flowcell, context={'request': request})

        data = serializer.data

        return_dict = {"data": data}

        return Response(return_dict)

    elif request.method == 'POST':

        serializer = FlowcellSerializer(data=request.data, context={'request': request})

        if serializer.is_valid():

            serializer.save()

            return Response(serializer.data, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET'])
def flowcell_summary_barcode(request, pk):

    q = request.GET.get('barcode_name', 'All reads')

    flowcell = Flowcell.objects.get(pk=pk)

    qs = FlowcellSummaryBarcode.objects \
        .filter(flowcell=flowcell) \
        .filter(barcode_name=q) \
        # .filter(run_id__owner=request.user)

    serializer = FlowcellSummaryBarcodeSerializer(qs, many=True)

    return Response(serializer.data)


def flowcell_summary_html(request, pk):

    flowcell = Flowcell.objects.get(pk=pk)

    qs = FlowcellSummaryBarcode.objects \
        .filter(flowcell=flowcell) \
        .exclude(barcode_name='No barcode')

    total = 0
    for row in qs:
        if row.barcode_name=="All reads":
            total+=row.total_length

    barcodedict=dict()

    for row in qs:
        if row.read_type_name not in barcodedict:
            barcodedict[row.read_type_name]=dict()
        if row.barcode_name not in barcodedict[row.read_type_name]:
            barcodedict[row.read_type_name][row.barcode_name]=dict()
        barcodedict[row.read_type_name][row.barcode_name][row.status]=1
        row.percentage = round(row.total_length/total*100,2)
        if row.read_count>0:
            row.avg_qual = round(row.quality_sum/row.read_count,2)
        else:
            row.avg_qual = 0

    """
    qs = list(qs)
    tempqs = deepcopy(qs[-1])
    for readtype in barcodedict:
        for barcode in barcodedict[readtype]:
            if len(barcodedict[readtype][barcode])<2:
                #print (type(barcodedict[readtype][barcode]))
                for status in barcodedict[readtype][barcode]:
                    if status==True:
                        tempqs.read_type_name=readtype
                        tempqs.barcode_name=barcode
                        tempqs.read_count=0
                        tempqs.total_length=0
                        print (tempqs.status)
                        tempqs.status=None
                        print (tempqs.status)
                        qs.append(tempqs)
                        pass
                    else:
                        tempqs.read_type_name=readtype
                        tempqs.barcode_name=barcode
                        tempqs.read_count = 0
                        tempqs.total_length = 0
                        print (tempqs.status)
                        tempqs.status=None
                        print (tempqs.status)
                        qs.append(tempqs)
                        pass

    #qs.append(qs[-1])
    print (qs)
    
    """

    ### Manipulate qs to add in missing records?!

    return render(request, 'reads/flowcell_summary.html', {'qs': qs})


@api_view(['GET'])
def flowcell_statistics(request, pk):
    """
    Return a prepared set of summaries for reads quality over time grouped by minute.
    """

    q = request.GET.get('barcode_name', 'All reads')

    flowcell = Flowcell.objects.get(pk=pk)

    run_list = flowcell.runs.all()

    key_list = FlowcellStatisticBarcode.objects\
        .filter(flowcell=flowcell)\
        .filter(barcode_name=q)\
        .values_list('barcode_name', 'read_type_name', 'status')\
        .distinct()

    result_dict = dict()

    data_keys = []

    for (barcode_name, read_type_name, status) in key_list:

        queryset = FlowcellStatisticBarcode.objects\
            .filter(flowcell=flowcell)\
            .filter(barcode_name=barcode_name)\
            .filter(read_type_name=read_type_name)\
            .filter(status=status)\
            .order_by('sample_time')

        if len(queryset) > 0:

            quality_list = []

            cumulative_bases = 0
            cumulative_reads = 0

            for record in queryset:

                if status == 'True':
                    l_status = 'Pass'

                else:
                    l_status = 'Fail'

                cumulative_bases = cumulative_bases + record.total_length
                cumulative_reads = cumulative_reads + record.read_count

                quality_list.append((
                    record.sample_time.timestamp() * 1000,
                    float(record.quality_sum/record.read_count),  # chart average quality over time
                    float(record.total_length / record.read_count),  # chart average read length over time
                    int(cumulative_bases),  # chart cumulative bases
                    int(cumulative_reads),  # chart cumulative reads
                    int(record.max_length),  # chart
                    float(record.total_length/60),  # chart sequence rate
                    float(record.total_length/record.channel_count/60),  # chart sequence speed
                ))

            result_dict["{} - {} - {}".format(barcode_name, read_type_name, l_status)] = quality_list

            data_keys.append("{} - {} - {}".format(barcode_name, read_type_name, l_status))

    run_data = []

    for run in run_list:
        run_dict = {
            'id': run.id,
            'name': run.name,
            'runid': run.runid,
            'start_time': run.start_time,
            'last_read': run.last_read(),
        }

        run_data.append(run_dict)

    res = {
        'runs': run_data,
        'data_keys': data_keys,
        'data': result_dict,
        'date_created': datetime.datetime.now()
    }

    return HttpResponse(json.dumps(res, cls=DjangoJSONEncoder), content_type="application/json")


@api_view(['GET'])
def flowcell_histogram_summary(request, pk):

    barcode_name = request.GET.get('barcode_name', 'All reads')

    flowcell = Flowcell.objects.get(pk=pk)

    max_bin_index = FlowcellHistogramSummary.objects \
        .filter(flowcell=flowcell) \
        .filter(barcode_name=barcode_name) \
        .aggregate(Max('bin_index'))

    key_list = FlowcellHistogramSummary.objects\
        .filter(flowcell=flowcell)\
        .filter(barcode_name=barcode_name)\
        .values_list('barcode_name', 'read_type_name', 'status')\
        .distinct()

    data_keys = {}

    for (barcode_name, read_type_name, status) in key_list:

        queryset = FlowcellHistogramSummary.objects\
            .filter(flowcell=flowcell)\
            .filter(barcode_name=barcode_name)\
            .filter(read_type_name=read_type_name)\
            .filter(status=status)\
            .order_by('bin_index')

        result_read_count_sum = [0] * (max_bin_index['bin_index__max'] + 1)
        result_read_length_sum = [0] * (max_bin_index['bin_index__max'] + 1)

        if len(queryset) > 0:

            if status == 'True':
                l_is_pass = 'Pass'
            else:
                l_is_pass = 'Fail'

            l_key = ("{} - {} - {}".format(barcode_name, read_type_name, l_is_pass))

            for r in queryset:

                l_bin_index = r.bin_index
                l_read_count_sum = r.read_count
                l_read_length_sum = r.read_length

                result_read_count_sum[l_bin_index] = l_read_count_sum
                result_read_length_sum[l_bin_index] = l_read_length_sum

            if barcode_name not in data_keys.keys():

                data_keys[barcode_name] = {}

            if read_type_name not in data_keys[barcode_name].keys():

                data_keys[barcode_name][read_type_name] = {}

            data_keys[barcode_name][read_type_name][l_is_pass] = {
                'read_count': {
                    'name': l_key,
                    'data': result_read_count_sum
                },
                'read_length': {
                    'name': l_key,
                    'data': result_read_length_sum
                }
            }

    categories = list(range(1 * 900, (max_bin_index['bin_index__max'] + 1) * 900, 900))

    return Response({'data': data_keys, 'categories': categories})


@api_view(['GET'])
def flowcell_channel_summary(request, pk):
    """
    This function generates the json data for the heatmap charts
    (Reads per channel and Bases per channel). The javascript
    file requestHistogramData.js calls this function.
    """

    flowcell = Flowcell.objects.get(pk=pk)

    qs = FlowcellChannelSummary.objects.filter(flowcell=flowcell)

    result_mapped_to_flowcell = {}

    #
    # set all channels with value 0
    #

    result_mapped_to_flowcell = {}

    for i in range(1, flowcell.size + 1):

        coordinate = get_coords(i, flowcell.size)

        result_mapped_to_flowcell[i] = ([coordinate[0], coordinate[1], 0, 0])

    #
    # fill the channels with correct value
    #

    for record in qs:

        coordinate = get_coords(record.channel, flowcell.size)

        result_mapped_to_flowcell[record.channel] = ([coordinate[0], coordinate[1], record.read_count, record.read_length])

    return HttpResponse(json.dumps(list(result_mapped_to_flowcell.values()), cls=DjangoJSONEncoder), content_type="application/json")


@api_view(['GET'])
def flowcell_run_summaries_html(request, pk):

    flowcell = Flowcell.objects.get(pk=pk)

    result = []

    for run in flowcell.runs.all():

        # Create an element with information from the run (originally from the fastq files)

        element = {

            'runid': None,
            'run_start_time': None,
            'first_read': None,
            'last_read': None,
            'minknow_computer_name': None,
            'minion_id': None,
            'asic_id': None,
            'sequencing_kit': None,
            'purpose': None,
            'minknow_version': None,
            'flowcell_type': None,
            'flowcell_id': None,
            'sample_name': None,
            'experiment_name': None
        }

        minion_run_status_list = MinIONRunStatus.objects.filter(run_id=run)

        if len(minion_run_status_list) > 0:

            minion_run_status = minion_run_status_list[0]

            element['runid'] = minion_run_status.run_id.runid
            element['minknow_computer_name'] = minion_run_status.minKNOW_computer
            element['minion_id'] = minion_run_status.minION.minION_name
            element['asic_id'] = minion_run_status.minKNOW_asic_id
            element['sequencing_kit'] = minion_run_status.sequencing_kit
            element['purpose'] = minion_run_status.minKNOW_exp_script_purpose
            element['minknow_version'] = minion_run_status.minKNOW_version
            element['flowcell_type'] = minion_run_status.flowcell_type
            element['flowcell_id'] = minion_run_status.minKNOW_flow_cell_id
            element['sample_name'] = minion_run_status.minKNOW_sample_name
            element['experiment_name'] = minion_run_status.experiment_id

        result.append(element)

    return render(request, 'reads/flowcell_runs_summary.html', {
        'run_list': result
    })


@api_view(['GET'])
def flowcell_run_basecalled_summary_html(request, pk):

    flowcell = Flowcell.objects.get(pk=pk)

    result_basecalled_summary = []

    for run in flowcell.runs.all():

        try:

            run_summary = {
                'runid': run.summary.run.runid,
                'read_count': run.summary.read_count,
                'total_read_length': run.summary.total_read_length,
                'max_read_length': run.summary.max_read_length,
                'min_read_length': run.summary.min_read_length,
                'avg_read_length': run.summary.avg_read_length,
                'first_read_start_time': run.summary.first_read_start_time,
                'last_read_start_time': run.summary.last_read_start_time
            }

            result_basecalled_summary.append(run_summary)

        except RunSummary.DoesNotExist:

            print('RunSummary does not exist for run {}'.format(run.id))

    return render(request, 'reads/flowcell_run_basecalled_summary.html', {
        'result_basecalled_summary': result_basecalled_summary
    })


@api_view(['GET'])
def flowcell_run_stats_latest(request, pk, checkid):
    """
    TODO describe function
    """
    flowcell = Flowcell.objects.get(pk=pk)

    minion_run_status_list = MinIONRunStatus.objects.filter(run_id__flowcell=flowcell)

    if minion_run_status_list.count() > 0:

        minion_run_status = minion_run_status_list[0]

    else:

        minion_run_status = None

    crazyminIONrunstats = MinIONRunStats.objects.filter(run_id__in=flowcell.runs.all(), id__gt=checkid)[:1000]

    result = []

    for minion_run_stats in crazyminIONrunstats:

        element = {
            "id": minion_run_stats.id,
            "minION": None,
            "run_id": None,
            "sample_time": minion_run_stats.sample_time,
            "event_yield": minion_run_stats.event_yield,
            "asic_temp": minion_run_stats.asic_temp,
            "heat_sink_temp": minion_run_stats.heat_sink_temp,
            "voltage_value": minion_run_stats.voltage_value,
            "mean_ratio": minion_run_stats.mean_ratio,
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
            "occupancy": minion_run_stats.occupancy(),
            "minKNOW_read_count": minion_run_stats.minKNOW_read_count,
            "minKNOW_histogram_values": minion_run_stats.minKNOW_histogram_values,
            "minKNOW_histogram_bin_width": minion_run_stats.minKNOW_histogram_bin_width
        }

        if minion_run_status:

            element['minKNOW_colours_string'] = minion_run_status.minKNOW_colours_string

        else:

            element['minKNOW_colours_string'] = None

        result.append(element)

    return JsonResponse(result, safe=False)


@api_view(['GET'])
def flowcell_tasks_detail_all(request, pk):

    flowcell = Flowcell.objects.get(pk=pk)

    run_list = Run.objects.filter(flowcell=flowcell)

    queryset = JobType.objects.filter(private=False)

    result = []

    for jobtype in queryset:

        obj = {}
        obj.update({
            'name': jobtype.name,
            'description': jobtype.description,
            'long_description': jobtype.long_description,
            'read_count': jobtype.readcount,
            'reference': jobtype.reference,
            'transcriptome': jobtype.transcriptome
        })

        jobmasterlist = JobMaster.objects.filter(Q(run__in=run_list) | Q(flowcell=flowcell)).filter(job_type=jobtype)

        if len(jobmasterlist) > 0:
            obj2 = {}
            if jobmasterlist[0].reference:
                reference_name = jobmasterlist[0].reference.reference_name
            else:
                reference_name = ''
            obj2.update({
                'reference': reference_name,
                'last_read': jobmasterlist[0].last_read,
                'read_count': jobmasterlist[0].read_count,
                'temp_file': jobmasterlist[0].tempfile_name,
                'complete': jobmasterlist[0].complete,
                'running': jobmasterlist[0].running
            })

            if len(obj2) > 0:

                obj.update({
                    'job_details': obj2
                })

        result.append(obj)

    return HttpResponse(json.dumps(result), content_type="application/json")


@api_view(['GET'])
def flowcell_tabs_details(request, pk):
    """
    Return tab_id, tab_title, and tab_position for a given flowcell.
    """
    flowcell_tabs_dict = {
        "LiveEvent": {
            "id": "tab-live-event-data",
            "title": "Live Event Data",
            "position": 1
        },
        "ChanCalc": {
            "id": "tab-basecalled-data",
            "title": "Basecalled Data",
            "position": 2
        },
        "Reads": {
            "id": "tab-reads",
            "title": "Reads data",
            "position": 3
        },
        "Minimap2": {
            "id": "tab-sequence-mapping",
            "title": "Sequence Mapping",
            "position": 4
        },
        "Assembly": {
            "id": "tab-sequence-assembly",
            "title": "Assembly",
            "position": 5
        },
        "Minimap2_trans": {
            "id": "tab-transcriptome-mapping",
            "title": "Transcriptome Mapping",
            "position": 6
        },
        "ExportReads": {
            "id": "tab-export-reads",
            "title": "Download Read Data",
            "position": 7
        },
        "Runs": {
            "id": "tab-runs",
            "title": "Runs",
            "position": 8
        },
        "Metagenomics": {
            "id": "tab-metagenomics",
            "title": "Metagenomics",
            "position": 9
        }
    }

    tabs = ['summary-data', 'tasks']

    #
    # Check for basecalled data
    #
    flowcell_summary_barcode_list = FlowcellSummaryBarcode.objects.filter(flowcell__id=pk)

    if flowcell_summary_barcode_list.count() > 0:

        tabs.append('basecalled-data')
        tabs.append('reads')

    minion_run_stats_list = MinIONRunStats.objects.filter(run_id__flowcell__id=pk)

    if minion_run_stats_list.count() > 0:

        tabs.append('live-event-data')

    paf_rough_cov_list = PafRoughCov.objects.filter(job_master__flowcell_id=pk)

    if paf_rough_cov_list.count() > 0:

        tabs.append('sequence-mapping')

    centrifuge_output_list = CentrifugeOutput.objects.filter(task__flowcell_id=pk)

    if centrifuge_output_list.count() > 0:

        tabs.append('metagenomics')

    # tabs_send = list()
    #
    # flowcell = Flowcell.objects.get(pk=pk)
    #
    # run_list = Run.objects.filter(flowcell=flowcell)
    #
    # if MinIONRunStatus.objects.filter(run_id__in=run_list):
    #
    #     tabs.append(flowcell_tabs_dict['LiveEvent'])
    #
    # for master in JobMaster.objects.filter(Q(run__in=run_list) | Q(flowcell=flowcell)).filter(last_read__gt=0).values_list('job_type__name', flat=True):
    #     if master in flowcell_tabs_dict.keys():
    #
    #         tabs.append(flowcell_tabs_dict[master])
    #
    #     else:
    #
    #         print("Flowcell '" + pk + "' has JobType '" + master + "' but there is no corresponding tab defined in reads/views.py")
    #
    # for tab in tabs:
    #
    #     if tab not in tabs_send:
    #
    #         tabs_send.append(tab)

    return Response(tabs)


@api_view(['GET', 'POST'])
def grouprun_list(request):

    if request.method == 'GET':

        queryset = GroupRun.objects.distinct().filter(owner=request.user)

        serializer = GroupRunSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)

    elif request.method == 'POST':

        serializer = GroupRunSerializer(data=request.data, context={'request': request})

        if serializer.is_valid():

            serializer.save(owner=request.user)

            return Response(serializer.data, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET'])
def grouprun_detail(request, pk):

    search_criteria = request.GET.get('search_criteria', 'id')

    if search_criteria == 'id':

        grouprun_list = GroupRun.objects.filter(owner=request.user).filter(id=pk)

    elif search_criteria == 'name':

        grouprun_list = GroupRun.objects.filter(owner=request.user).filter(name=pk)

    else:

        grouprun_list = GroupRun.objects.none()

    if len(grouprun_list) != 1:

        return Response(status=status.HTTP_404_NOT_FOUND)

    grouprun = grouprun_list[0]

    serializer = GroupRunSerializer(grouprun, context={'request': request})

    return Response(serializer.data)


@api_view(['GET', 'POST'])
def grouprun_membership_list(request):

    if request.method == 'GET':

        grouprun_id = request.GET.get('grouprun_id')

        queryset = GroupRun.objects.filter(id=grouprun_id).filter(owner=request.user)

        if len(queryset) == 1:
            serializer = RunSerializer(queryset.first().runs, many=True, context={'request': request})

            return Response(serializer.data)

    elif request.method == 'POST':

        grouprun_id = request.data['grouprun_id']
        run_id = request.data['run_id']

        grouprun = GroupRun.objects.get(pk=grouprun_id)
        run = Run.objects.get(pk=run_id)

        run.groupruns.add(grouprun)
        run.save()

        serializer = RunSerializer(run, context={'request': request})

        return Response(serializer.data, status=status.HTTP_201_CREATED)


@api_view(['GET'])
def read_detail(request,pk):
    return Response({})


@api_view(['GET', 'POST'])
def read_list_new(request):

    search_criteria = request.GET.get('search_criteria', 'name')

    search_value = request.GET.get('search_value', 'name')

    offset = int(request.GET.get('offset', '0'))

    limit = int(request.GET.get('limit', '10'))

    if request.method == 'GET':

        if search_criteria == 'run':

            qs = FastqRead.objects.filter(run__id=search_value)[offset:offset + limit]

        elif search_criteria == 'fastqfile':

            qs = FastqRead.objects.filter(fastqfile_id=search_value)[offset:offset + limit]


        else:

            qs = FastqRead.objects.none()

        serializer = FastqReadGetSerializer(qs, many=True, context={'request': request})

        return Response(serializer.data)

    elif request.method == 'POST':

        # serializer = FastqReadSerializer(data=request.data, many=True)

        # print(request.data)

        logger.info('>>> received reads post - calling task - request.data size: {}'.format(len(request.data)))

        # save_reads.delay(request.data)
        save_reads(request.data)

        return Response({}, status=status.HTTP_201_CREATED)

        # if serializer.is_valid():
        #
        #     print('reads valid')
        #     res = serializer.save()
        #     print(res)
        #
        #     return Response({}, status=status.HTTP_201_CREATED)
        #
        # else:
        #
        #     print('reads not valid')
        #
        # return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

@task
def save_reads(data):

    serializer = FastqReadSerializer(data=data, many=True)

    if serializer.is_valid():
        serializer.save()
        logger.info('Saving reads with success')

    else:
        logger.info('Saving reads with failure')


@api_view(['GET'])
def readextra_list(request):

    search_criteria = request.GET.get('search_criteria', 'read_id')

    search_value = request.GET.get('search_value', 'name')

    if search_criteria == 'read_id':

        qs = FastqRead.objects.filter(read_id=search_value)

        if len(qs) > 0:

            fastq_record = qs[0]

            fastq_string = ">{} runid={} ch={} read={} start_time={}\n{}\n".format(
                fastq_record.read_id,
                fastq_record.run.runid,
                fastq_record.channel,
                fastq_record.read,
                fastq_record.start_time.strftime("%Y-%m-%dT%H:%M:%SZ"),
                fastq_record.fastqreadextra.sequence
            )

            read = {
                'sequence': fastq_string
            }

        else:
            read = {}

        return JsonResponse(read)


@api_view(['POST', 'GET'])
def barcode_list_new(request):

    if request.method == 'GET':

        barcode_list = Barcode.objects.all()

        serializer = BarcodeSerializer(barcode_list, many=True, context={'request': request})

        return Response(serializer.data)

    elif request.method == 'POST':

        serializer = BarcodeSerializer(data=request.data, context={'request': request})

        if serializer.is_valid():

            serializer.save()

            return Response(serializer.data, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


def version(request):

    version = {
        'server': '1.0',
        'clients': [
            '1.0'
        ],
        'minknow': [
            '1.11.5',
        ]
    }

    return HttpResponse(json.dumps(version), content_type="application/json")


class FlowcellTabList(APIView):

    def get(self, request, pk):

        flowcell_tab_list = FlowcellTab.objects.all()

        data = FlowcellTabSerializer(flowcell_tab_list, many=True).data

        return Response(data)
