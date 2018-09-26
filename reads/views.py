import datetime
import json
from datetime import timedelta

import dateutil.parser
import math
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

from jobs.models import JobMaster, JobType
from minotourapp import settings
from reads.models import (Barcode, FastqRead, FastqReadType,
                          MinIONControl, MinIONEvent,
                          MinIONEventType, MinionMessage, MinIONRunStats,
                          MinIONRunStatus, MinIONScripts, MinIONStatus, Run, GroupRun, FlowcellStatisticBarcode,
                          FlowcellSummaryBarcode, Flowcell, MinION)
from reads.models import FlowcellChannelSummary
from reads.models import FlowcellHistogramSummary
from reads.serializers import (BarcodeSerializer,
                               ChannelSummarySerializer, FastqReadSerializer,
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
                               GroupRunSerializer, FlowcellSummaryBarcodeSerializer)
from reads.utils import get_coords
from reference.models import ReferenceInfo


def humanbases(n):
    """
    :purpose: Convert bases to human readable format
    :used_by: used by many
    :author: Matt Loose

    ChangeLog
    2018-07-09 Add documentation

    :param n: (integer) number of bases
    :return: (string) human readable format
    """
    #return (n)
    millnames = ['', ' Kb', ' Mb', ' Gb', ' Tb']
    n = float(n)
    millidx = max(0,min(len(millnames)-1,
                        int(math.floor(0 if n == 0 else math.log10(abs(n))/3))))

    return '{:.3f}{}'.format(n / 10**(3 * millidx), millnames[millidx])




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
    if request.method == 'GET':
        queryset = FastqReadType.objects.all()
        serializer = FastqReadTypeSerializer(queryset, many=True, context={'request': request})
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
    if request.method == "GET":
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

    if request.method == "GET":
        serializer = MinIONEventTypeSerializer(event_, context={'request': request})
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

    if request.method == 'GET':
        serializer = FastqReadTypeSerializer(run, context={'request': request})
        return Response(serializer.data)


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
        #print ("run list ", request.user)
        queryset = Run.objects.filter(owner=request.user).filter(to_delete=False)
        serializer = RunSerializer(queryset, many=True, context={'request': request})
        #for q in queryset:
        #    print ("run queryset",q.owner)
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



@api_view(['GET'])
def current_run_list(request):
    """
    :purpose: Get returns a list of all runs in minotour owned by a specific user which are specified as active.
    :used_by: minotour app uses this endpoint to identify active runs.
    :author: Roberto Santos

    ChangeLog
    2018-07-09 Add documentation - Matt

    :param request: (standard django request) without querystring parameter
    :return: (str) json format
    """
    if request.method == 'GET':
        #queryset = MinIONRun.objects.filter(owner=request.user).filter(Q(reads__created_date__gte = datetime.now()-timedelta(days=1))  | Q(RunStats__created_date__gte = datetime.now()-timedelta(days=1) )).distinct()
        #print ("current_run:",request.user)
        queryset = Run.objects.filter(owner=request.user).filter(active=True).distinct()
        #for q in queryset:
        #    print ("current queryset:", q.owner)
        serializer = RunSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)


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
    if request.method == 'GET':
        #queryset = MinION.objects.filter(owner=request.user)
        queryset = [obj for obj in MinION.objects.filter(owner=request.user) if obj.status() != "unplugged"]
        serializer = MinIONSerializer(queryset, many=True, context={'request': request})
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
def sinceminion_messages_list(request, pk, starttime,endtime):
    if request.method == 'GET':
        #print (starttime)
        correctedstart = parser.parse(starttime) - timedelta(minutes=180)
        correctedend = parser.parse(endtime) + timedelta(minutes=180)
        #print (correctedstart.isoformat().replace('+00:00', 'Z'))
        queryset = MinionMessage.objects.filter(minION=pk).filter(
        minKNOW_message_timestamp__gte=correctedstart.isoformat().replace('+00:00', 'Z')).filter(minKNOW_message_timestamp__lte=correctedend.isoformat().replace('+00:00', 'Z'))
        serializer = MinionMessageSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)


@api_view(['GET'],)
def recentminion_messages_list(request, pk):

    if request.method == 'GET':
        queryset = MinionMessage.objects.filter(minION=pk).filter(minKNOW_message_timestamp__gte=timezone.now() - timedelta(hours=24))
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

    if form_start_time:
        start_time = dateutil.parser.parse(form_start_time)

    else:
        start_time = None

    if form_end_time:
        end_time = dateutil.parser.parse(form_end_time)

    else:
        end_time = None

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
            .filter(minion__in=minion_list)[:10] \
            .order_by('-timestamp')

    serializer = MinionMessageSerializer(messages, many=True, context={'request': request})

    return Response(serializer.data)


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
        #print (serializer)
        if serializer.is_valid():
            serializer.save(owner=request.user)
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET','POST'],)
def minION_control_update(request,pk,checkid):
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
        #print (serializer)
        return Response(serializer.data)

    if request.method == 'POST':
        #serializer = MinIONControlSerializer(data=event_, context={'request': request})
        #if serializer.is_valid():
        #    serializer.complete=True
        #    serializer.save(update_fields=['complete'])

        event_ = MinIONControl.objects.get(id=checkid)
        event_.complete=True
        event_.save(update_fields=['complete'])
        return Response(status=status.HTTP_204_NO_CONTENT)
        #return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
        #
        # serializer = MinIONControlSerializer(event_, context={'request': request})
        #if serializer.is_valid():
        #    serializer.complete=True
        #    serializer.save()
        #    return Response(serializer.data, status=status.HTTP_201_CREATED)
        #return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    #try:
    #    event_ = MinIONControl.objects.get(id=checkid)
    #except MinIONControl.DoesNotExist:
    #    return Response(status=status.HTTP_404_NOT_FOUND)


@api_view(['GET',])
def minION_currentrun_list(request, pk):
    """
    TODO describe function
    """
    try:
        #print (pk)
        #print (MinIONEvent.objects.all())
        minion = MinION.objects.get(pk=pk)
    except MinION.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        serializer = MinIONSerializer(minion, context={'request': request})
        return Response(serializer.data)


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
        #minIONrunstats = MinIONRunStats.objects.all()
        #print (len(crazyminIONrunstats))

    except MinIONRunStats.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        #print (crazyminIONrunstats)
        #for item in crazyminIONrunstats:
        #    print (item.minION, item.run_id, item.sample_time, item.event_yield)

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
        #minIONrunstats = MinIONRunStats.objects.all()
        #print (crazyminIONrunstats)

    except MinIONRunStats.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        #print (crazyminIONrunstats)
        #for item in crazyminIONrunstats:
        #    print (item.minION, item.run_id, item.sample_time, item.event_yield)

        serializer = MinIONRunStatsSerializer(crazyminIONrunstats, many=True , context={'request': request})

        return Response(serializer.data)


@api_view(['GET', 'PUT', 'POST', 'DELETE'])
def minION_run_status_list(request,pk):
    """
    TODO describe function
    """
    if request.method == 'POST':
        #print ("we are posting here")
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

        serializer = FastqReadSerializer(querysets, many=True, context={'request': request})

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
    # try:
    #     # print (pk)
    #     # print (MinIONEvent.objects.all())
    #     minion = MinION.objects.get(pk=pk)
    # except MinION.DoesNotExist:
    #     return Response(status=status.HTTP_404_NOT_FOUND)
    #
    # if request.method == 'GET':
    #     serializer = MinIONSerializer(minion, context={'request': request})
    #     return Response(serializer.data)

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

        print ("We are inside POST")
        print (request.data)

        # id of minION to update
        # parameter to update
        # new value

        minion = MinION.objects.get(pk=pk)

        serializer = MinIONSerializer(minion, data=request.data, partial=True,context={'request': request})

        if serializer.is_valid():
            serializer.save()

        else:
            print ("PROBLEM")

        return Response(serializer.data)


@api_view(['GET', 'PUT', 'DELETE'])
def read_detail(request, pk):
    """
    Retrieve, update or delete a read instance.
    """
    try:
        read = FastqRead.objects.get(pk=pk)
    except FastqRead.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        serializer = FastqReadSerializer(read, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'PUT':
        serializer = FastqReadSerializer(read, data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    elif request.method == 'DELETE':
        read.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


"""
@api_view(['PATCH'])
def read_update(request,pk,readid):
    if request.method == "PATCH":
        print (request.data)
        print (readid)
        FastQRecord = FastqRead.objects.get(read_id=readid)
        #FastQRecord = FastqRead.objects.filter(run_id=pk)
        print (FastQRecord.type_id)
        print (request.data["type"])
        if (str(FastQRecord.type_id) != request.data["type"].split('/')[-2]):
            serializer = FastqReadSerializer(FastQRecord, data=request.data, context={'request': request}, partial=True)
            if serializer.is_valid():
                serializer.save()
                FastQRecord = FastqRead.objects.get(read_id=readid)
                try:
                    update_global_state(FastQRecord)
                except Exception as e:
                    print ("update global state failed", e)
                return Response(serializer.data, status=status.HTTP_201_CREATED)
            print('Serializer errors: {} '.format(serializer.errors))
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
        else:
            print ('State not needed to change.')
            return Response("state not needed", status=status.HTTP_201_CREATED)
"""


@api_view(['GET'])
def readname_list(request, pk):
    """
    TODO describe function
    CUrrent API behaviour asssumes starting from page1.
    Request for page 0 returns last page.
    This could be a source of confusion and we should resolve.
    """
    if request.method == 'GET':
        queryset = FastqRead.objects.filter(run_id=pk).order_by('id')

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
def cumulative_read_count(request, pk):
    """
    TODO describe function
    """
    if request.method == 'GET':
        queryset = RunStatistic.objects.filter(run_id__owner=request.user).filter(run_id=pk).order_by('type',
                                                                                                      'sample_time', )
        result = dict()
        for key in queryset:
            if str(key.type) not in result:
                result[str(key.type)] = dict()
            result[str(key.type)][str(key.sample_time)] = key.total_length
        data_to_return = dict()
        counter = 0
        for readtype in result:
            cumuyield = 0
            data_to_return[counter] = dict()
            data_to_return[counter]['name'] = readtype
            data_to_return[counter]['data'] = list()
            for data in result[readtype]:
                cumuyield = cumuyield + result[readtype][data]
                data_to_return[counter]['data'].append((UTC_time_to_epoch(data), cumuyield))
        return HttpResponse(json.dumps(data_to_return), content_type="application/json")


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
        print (request)

        #run=request.POST["run"]
        barcodename=request.data["name"]
        print (request.data)
        print(barcodename)

        minionrun = Run.objects.get(pk=pk)

        flowcellruns = minionrun.flowcellrun

        print (flowcellruns)

        barcode,created2 = Barcode.objects.get_or_create(run=minionrun, name=barcodename)

        serializer = BarcodeSerializer(barcode, context={'request': request})

        #if serializer.is_valid():
        #if created2:
        return Response(serializer.data, status=status.HTTP_201_CREATED)

        #return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET'])
def barcode_detail(request, pk):
    """
    Return the details of a particular barcode.
    """
    queryset = Barcode.objects.get(pk=pk)

    serializer = BarcodeSerializer(queryset, many=False, context={'request': request})

    return Response(serializer.data)


def UTC_time_to_epoch(timestamp):
    """
    TODO describe function
    """
    dt = parser.parse(timestamp)

    return dt.timestamp()*1000


@api_view(['GET'])
def minION_liverun_list(request,pk):
    """
    TODO describe function
    """
    return None


@api_view(['GET'])
def flowcell_list_active(request):

    queryset = Flowcell.objects.distinct().filter(owner=request.user)

    serializer = FlowcellSerializer(queryset, many=True, context={'request': request})

    return Response(serializer.data)


@api_view(['GET','POST'])
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
                'number_runs': record.number_runs,
                'number_barcodes': record.number_barcodes,
                'total_read_length': record.total_read_length,
                'average_read_length': record.average_read_length,
                'is_active': record.is_active,
                'sample_name': record.sample_name
            }

            flowcells.append(flowcell)

        # serializer = FlowcellSerializer(queryset, many=True, context={'request': request})
        # return Response(serializer.data)
        return JsonResponse({'data': flowcells})

    elif request.method == 'POST':
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

            return Response(status=status.HTTP_404_NOT_FOUND)

        flowcell = flowcell_list[0]

        serializer = FlowcellSerializer(flowcell, context={'request': request})

        return Response(serializer.data)

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
def flowcell_run_status_list(request, pk):

    flowcell = Flowcell.objects.get(pk=pk)
    queryset = MinIONRunStatus.objects.filter(run_id__in=flowcell.runs.all())
    #print (queryset)
    serializer = MinIONRunStatusSerializer(queryset, many=True, context={'request': request})
    return Response(serializer.data)


@api_view(['GET'])
def flowcell_run_stats_latest(request, pk, checkid):
    """
    TODO describe function
    """
    flowcell = Flowcell.objects.get(pk=pk)
    # queryset = FlowCellRun.objects.filter(flowcell_id=pk)
    # runset = list()
    # for run in queryset:
    #     # print (run.run_id)
    #     runset.append(run.run_id)

    try:
        crazyminIONrunstats = MinIONRunStats.objects.filter(run_id__in=flowcell.runs.all(), id__gt=checkid)[:1000]
        #minIONrunstats = MinIONRunStats.objects.all()
        #print (crazyminIONrunstats)

    except MinIONRunStats.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        #print (crazyminIONrunstats)
        #for item in crazyminIONrunstats:
        #    print (item.minION, item.run_id, item.sample_time, item.event_yield)

        serializer = MinIONRunStatsSerializer(crazyminIONrunstats, many=True , context={'request': request})

        return Response(serializer.data)


@api_view(['POST'])
def flowcellset_task_detail_all(request,pk):
    """We need to check if a job type already exists - if so we are not going to add another."""
    if request.method == 'POST':
        jobtype=JobType.objects.get(name=request.data["job"])
        print (jobtype)
        reference=""
        if request.data["reference"]!="null":
            reference=ReferenceInfo.objects.get(reference_name=request.data["reference"])
            print (reference)
        flowcellitem=Flowcell.objects.get(id=pk)
        print (flowcellitem)
        print(request.data)
        print(jobtype,reference,flowcellitem)
        jobmasters=JobMaster.objects.filter(flowcell=flowcellitem).filter(job_type=jobtype)
        print ("Jobmasters",jobmasters)
        if len(jobmasters) > 0:
            #return Response("Duplicate Job attempted. Not allowed.", status=status.HTTP_400_BAD_REQUEST)
            return Response("Duplicate Job attempted. Not allowed.", status=status.HTTP_200_OK)
        else:
            newjob = JobMaster(flowcell=flowcellitem, job_type=jobtype, last_read=0, read_count=0, complete=False, running=False)

            print ("trying to make a job", newjob)

            if request.data["reference"] != "null":
            #if len(reference)>0:
                newjob.reference = reference
            try:
                newjob.save()
                print ("job created")
            except Exception as e:
                print (e)
                print ("error")
            return Response("Job Created.", status=status.HTTP_200_OK)


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
def tabs_details(request, pk):
    """
    Return tab_id, tab_title, and tab_position for a given run.
    """
    run_tabs_dict = {
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
        "Metagenomics": {
            "id": "tab-metagenomics",
            "title": "Sequence Identification",
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
        }
    }
    tabs = list()
    # Find live event data
    if MinIONRunStatus.objects.filter(run_id=pk):
        tabs.append(run_tabs_dict['LiveEvent'])

    for master in JobMaster.objects.filter(run_id=pk).values_list('job_type__name', flat=True):
        if master in run_tabs_dict.keys():
            tabs.append(run_tabs_dict[master])
        else:
            print("RunID '" + pk + "' has JobType '" + master + "' but there is no corresponding tab defined in reads/views.py")

    return Response(tabs)


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
        "Kraken": {
            "id": "tab-sequence-id",
            "title": "Sequence Identification",
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

    tabs = list()

    tabs_send = list()

    flowcell = Flowcell.objects.get(pk=pk)
    run_list = Run.objects.filter(flowcell=flowcell)

    if MinIONRunStatus.objects.filter(run_id__in=run_list):

        tabs.append(flowcell_tabs_dict['LiveEvent'])

    for master in JobMaster.objects.filter(Q(run__in=run_list) | Q(flowcell=flowcell)).filter(last_read__gt=0).values_list('job_type__name', flat=True):

        if master in flowcell_tabs_dict.keys():

            tabs.append(flowcell_tabs_dict[master])

        else:

            print("Flowcell '" + pk + "' has JobType '" + master + "' but there is no corresponding tab defined in reads/views.py")

    for tab in tabs:

        if tab not in tabs_send:

            tabs_send.append(tab)

    # tabs_send.append(flowcell_tabs_dict['Runs'])  # always add a tab for runs
    # tabs_send.append(flowcell_tabs_dict['Metagenomics'])  # always add a tab for runs

    return Response(tabs_send)


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


@api_view(['GET', 'POST'])
def read_list_new(request):

    search_criteria = request.GET.get('search_criteria', 'name')

    search_value = request.GET.get('search_value', 'name')

    if request.method == 'GET':

        if search_criteria != 'run':

            qs = FastqRead.objects.filter(run__id=search_value)

        elif search_criteria == 'grouprun':

            qs = FastqRead.objects.filter(run__id=search_value)

        serializer = FastqReadSerializer(qs, many=True, context={'request': request})

        return Response(serializer.data)

    elif request.method == 'POST':

        serializer = FastqReadSerializer(data=request.data, many=True, context={'request': request})

        if serializer.is_valid():

            serializer.save()

            return Response({}, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET'])
def readextra_list(request):

    search_criteria = request.GET.get('search_criteria', 'read_id')

    search_value = request.GET.get('search_value', 'name')

    if search_criteria == 'read_id':

        qs = FastqRead.objects.filter(read_id=search_value)

        if len(qs) > 0:
            read = {
                'sequence': qs[0].fastqreadextra.sequence
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
