import json
from datetime import datetime, timedelta
from django.utils import timezone
from dateutil import parser
from django.db.models import Q
from django.http import HttpResponse
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response
#from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.core.serializers.json import DjangoJSONEncoder

import numpy as np

from reads.models import Barcode, ChannelSummary, BarcodeGroup
from reads.models import FastqRead
from reads.models import FastqReadExtra
from reads.models import FastqReadType
from reads.models import HistogramSummary
from reads.models import MinION
from reads.models import MinIONControl
from reads.models import MinIONEvent
from reads.models import MinIONEventType
from reads.models import MinIONRun
from reads.models import MinIONRunStats
from reads.models import MinIONRunStatus
from reads.models import MinIONScripts
from reads.models import MinIONStatus
from reads.models import MinIONmessages
from reads.models import RunStatistic
from reads.models import RunStatisticBarcode
from reads.models import RunSummary
from reads.models import RunSummaryBarcode
from reads.models import JobType
from reads.models import JobMaster
from reads.models import FlowCell
from reads.models import FlowCellRun
from reference.models import ReferenceInfo
from reads.serializers import BarcodeSerializer, BarcodeGroupSerializer
from reads.serializers import ChannelSummarySerializer
from reads.serializers import FastqReadSerializer
from reads.serializers import RunSummaryBarcodeSerializer
from reads.serializers import RunStatisticBarcodeSerializer
from reads.serializers import RunStatisticSerializer
from reads.serializers import FastqReadTypeSerializer
from reads.serializers import MinIONControlSerializer
from reads.serializers import MinIONEventSerializer
from reads.serializers import MinIONEventTypeSerializer
from reads.serializers import MinIONRunSerializer
from reads.serializers import MinIONRunStatsSerializer
from reads.serializers import MinIONRunStatusSerializer
from reads.serializers import MinIONScriptsSerializer
from reads.serializers import MinIONSerializer
from reads.serializers import MinIONStatusSerializer
from reads.serializers import MinIONmessagesSerializer
from reads.serializers import RunHistogramSummarySerializer
from reads.serializers import RunSummarySerializer
from reads.serializers import JobTypeSerializer
from reads.serializers import FlowCellSerializer
from reads.serializers import FlowCellRunSerializer

from reads.models import update_global_state,update_sum_stats

from minotourapp import settings

@api_view(['GET'])
def read_type_list(request):
    """
    List of all runs by fastqread types
    """
    if request.method == 'GET':
        queryset = FastqReadType.objects.all()
        serializer = FastqReadTypeSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)


@api_view(['GET'])
def events_type_list(request):
    """
    Lis of all event types.
    """
    if request.method == "GET":
        queryset = MinIONEventType.objects.all()
        serializer = MinIONEventTypeSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)


@api_view(['GET'])
def events_type_detail(request, pk):
    """
    Lis of all event types.
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
    Retrieve a fastqread type instance.
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
    List of all runs by user, or create a new run.
    """
    if request.method == 'GET':
        #print ("run list ", request.user)
        queryset = MinIONRun.objects.filter(owner=request.user).filter(to_delete=False)
        serializer = MinIONRunSerializer(queryset, many=True, context={'request': request})
        #for q in queryset:
        #    print ("run queryset",q.owner)
        return Response(serializer.data)

    elif request.method == 'POST':
        serializer = MinIONRunSerializer(data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save(owner=request.user)
            return Response(serializer.data, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

@api_view(['GET'])
def current_run_list(request):
    """
    List of all runs by user, or create a new run.
    """
    if request.method == 'GET':
        #queryset = MinIONRun.objects.filter(owner=request.user).filter(Q(reads__created_date__gte = datetime.now()-timedelta(days=1))  | Q(RunStats__created_date__gte = datetime.now()-timedelta(days=1) )).distinct()
        #print ("current_run:",request.user)
        queryset = MinIONRun.objects.filter(owner=request.user).filter(active=True).distinct()
        #for q in queryset:
        #    print ("current queryset:", q.owner)
        serializer = MinIONRunSerializer(queryset, many=True, context={'request': request})
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
def minION_messages_list(request, pk):
    """
    TODO describe function
    """
    if request.method == 'GET':
        queryset = MinIONmessages.objects.filter(minION=pk)
        serializer = MinIONmessagesSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'POST':
        serializer = MinIONmessagesSerializer(data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

@api_view(['GET'],)
def sinceminION_messages_list(request, pk,starttime,endtime):
    if request.method == 'GET':
        #print (starttime)
        correctedstart = parser.parse(starttime) - timedelta(minutes=180)
        correctedend = parser.parse(endtime) + timedelta(minutes=180)
        #print (correctedstart.isoformat().replace('+00:00', 'Z'))
        queryset = MinIONmessages.objects.filter(minION=pk).filter(
        minKNOW_message_timestamp__gte=correctedstart.isoformat().replace('+00:00', 'Z')).filter(minKNOW_message_timestamp__lte=correctedend.isoformat().replace('+00:00', 'Z'))
        serializer = MinIONmessagesSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)

@api_view(['GET'],)
def recentminION_messages_list(request, pk):
    """
    TODO describe function
    """
    if request.method == 'GET':
        queryset = MinIONmessages.objects.filter(minION=pk).filter(minKNOW_message_timestamp__gte= timezone.now() - timedelta(hours=24))
        serializer = MinIONmessagesSerializer(queryset, many=True, context={'request': request})
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
    try:
        run = MinIONRun.objects.get(pk=pk)

        if run.owner != request.user:
            return Response(status=status.HTTP_401_UNAUTHORIZED)

    except MinIONRun.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        serializer = MinIONRunSerializer(run, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'PUT':
        serializer = MinIONRunSerializer(run, data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    elif request.method == 'DELETE':
        run.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


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
        serializer = MinIONScriptsSerializer(data=request.data, context={'request': request})
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
    """
    List of all runs by user, or create a new run.
    """
    if request.method == 'GET':
        querysets = FastqRead.objects.filter(run_id=pk)
        #paginator = Paginator(querysets, 25)  # Show 25 contacts per page
        #page = page
        #print (paginator.count,paginator.num_pages)
        #try:
        #    queryset = paginator.page(page)
        #except PageNotAnInteger:
            # If page is not an integer, deliver first page.
        #    queryset = paginator.page(1)
        #except EmptyPage:
            # If page is out of range (e.g. 9999), deliver last page of results.
        #    queryset = paginator.page(paginator.num_pages)

        serializer = FastqReadSerializer(querysets, many=True, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'POST':
        serializer = FastqReadSerializer(data=request.data, many=True, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        #return Response('', status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET'])
def minION_detail(request, pk):
    """
    TODO describe function
    """
    try:
        # print (pk)
        # print (MinIONEvent.objects.all())
        minion = MinION.objects.get(pk=pk)
    except MinION.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        serializer = MinIONSerializer(minion, context={'request': request})
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
def run_summary(request, pk):
    """
    Return a list with summaries for each read type of a given run.
    """
    queryset = RunSummary.objects\
        .filter(run_id__owner=request.user)\
        .filter(run_id=pk)

    serializer = RunSummarySerializer(queryset, many=True, context={'request': request})

    return Response(serializer.data)


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
def run_summary_by_minute(request, pk, last=''):
    """
    Return a list with summaries for a particular run grouped by minute.
    """

    if last == '':
        queryset = RunStatistic.objects\
            .filter(run_id__owner=request.user)\
            .filter(run_id=pk)

    else:
        queryset = RunStatistic.objects\
            .filter(run_id__owner=request.user)\
            .filter(run_id=pk)\
            .filter(sample_time__gt=last)

    serializer = RunStatisticSerializer(queryset, many=True, context={'request': request})

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
def run_histogram_summary(request, pk):
    """
    Return a list with histogram summaries for a particular run.
    """
    queryset = HistogramSummary.objects\
        .filter(run_id__owner=request.user)\
        .filter(run_id=pk)

    serializer = RunHistogramSummarySerializer(queryset, many=True, context={'request': request})

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

        minionrun = MinIONRun.objects.get(pk=pk)

        flowcellruns = minionrun.flowcellrun

        print (flowcellruns)



        #Can we great the barcodegroup instance here?

        barcodegroup,created = BarcodeGroup.objects.get_or_create(flowcell=minionrun.flowcellrun.first().flowcell, name=barcodename)

        barcode,created2 = Barcode.objects.get_or_create(run=minionrun,barcodegroup=barcodegroup,name=barcodename)

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


@api_view(['GET'])
def barcodegroup_detail(request, pk):
    """
    Return the details of a particular barcode.
    """
    queryset = BarcodeGroup.objects.get(pk=pk)

    serializer = BarcodeGroupSerializer(queryset, many=False, context={'request': request})

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
def tasks_detail(request):
    if request.method == 'GET':
        queryset=JobType.objects.filter(private=False)
        serializer = JobTypeSerializer(queryset,many=True, context={'request': request})
        return Response(serializer.data)


@api_view(['POST'])
def set_task_detail_all(request,pk):
    """We need to check if a job type already exists - if so we are not going to add another."""
    if request.method == 'POST':
        jobtype=JobType.objects.get(name=request.data["job"])
        print (jobtype)
        reference=""
        if request.data["reference"]!="null":
            reference=ReferenceInfo.objects.get(reference_name=request.data["reference"])
            print (reference)
        minionrun=MinIONRun.objects.get(id=pk)
        print (minionrun)
        print(request.data)
        print(jobtype,reference,minionrun)
        jobmasters=JobMaster.objects.filter(run=minionrun).filter(job_type=jobtype)
        print ("Jobmasters",jobmasters)
        if len(jobmasters) > 0:
            #return Response("Duplicate Job attempted. Not allowed.", status=status.HTTP_400_BAD_REQUEST)
            return Response("Duplicate Job attempted. Not allowed.", status=status.HTTP_200_OK)
        else:
            newjob = JobMaster(run=minionrun, job_type=jobtype, last_read=0, read_count=0, complete=False, running=False)

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
def tasks_detail_all(request,pk):
    queryset=JobType.objects.filter(private=False).filter(type_name__type_name__in=['run',])
    minionrun=MinIONRun.objects.get(pk=pk)

    result = []
    #print (queryset)
    for jobtype in queryset:
        #print (jobtype.type_name.all())
        obj = {}
        obj.update({
            'name': jobtype.name,
            'description' : jobtype.description,
            'long_description' : jobtype.long_description,
            'read_count': jobtype.readcount,
            'reference': jobtype.reference,
            'transcriptome': jobtype.transcriptome
        })

        jobmasterlist = JobMaster.objects.filter(run=minionrun).filter(job_type=jobtype)

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

            obj.update({
                'job_details': obj2
            })

        result.append(obj)

    return HttpResponse(json.dumps(result), content_type="application/json")

@api_view(['GET'])
def flowcell_list_active(request):
    if request.method == 'GET':
        queryset = FlowCell.objects.distinct().filter(owner=request.user).filter(flowcelldetails__run__active=True)
        serializer = FlowCellSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)



@api_view(['GET','POST'])
def flowcell_list(request):
    #if request.method == 'GET':
    #    queryset = FlowCell.objects.filter(owner=request.user)
    #    serializer = FlowCellSerializer(queryset, many=True, context={'request': request})
    #    return Response(serializer.data)
    if request.method == 'GET':
        queryset = FlowCell.objects.distinct().filter(owner=request.user)
        serializer = FlowCellSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'POST':
        serializer = FlowCellSerializer(data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save(owner=request.user)
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET','POST'])
def flowcell_detail(request,pk):
    if request.method == 'GET':
        queryset = FlowCellRun.objects.filter(flowcell_id=pk)
        serializer = FlowCellRunSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)
    elif request.method == 'POST':
        serializer = FlowCellRunSerializer(data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

@api_view(['GET'])
def flowcell_summary_barcode(request,pk):
    queryset = FlowCellRun.objects.filter(flowcell_id=pk)
    runset = list()
    for run in queryset:
        #print (run.run_id)
        runset.append(run.run_id)
    queryset = RunSummaryBarcode.objects \
        .filter(run_id__owner=request.user) \
        .filter(run_id__in=runset)

    serializer = RunSummaryBarcodeSerializer(queryset, many=True, context={'request': request})

    return Response(serializer.data)

@api_view(['GET'])
def flowcell_summary_barcode_by_minute(request, pk):
    """
        Return a list with summaries for a particular run grouped by minute.
        """
    queryset = FlowCellRun.objects.filter(flowcell_id=pk)
    runset = list()
    for run in queryset:
        # print (run.run_id)
        runset.append(run.run_id)
    queryset = RunStatisticBarcode.objects \
        .filter(run_id__owner=request.user) \
        .filter(run_id__in=runset)

    serializer = RunStatisticBarcodeSerializer(queryset, many=True, context={'request': request})

    return Response(serializer.data)


@api_view(['GET'])
def flowcell_summary_barcode_by_minute_quality(request,pk):
    """
    Return a prepared set of summaries for reads quality over time grouped by minute.
    :param request:
    :param pk:
    :return:
    """
    barcodeset = BarcodeGroup.objects.filter(flowcell=pk)
    barcodedict = dict()
    resultset = dict()
    readtypes = FastqReadType.objects.all()
    for barcode in barcodeset:
        barcodedict[barcode.name] = barcode.id
        for readtype in readtypes:
            queryset = RunStatisticBarcode.objects.filter(barcode__barcodegroup=barcode.id).filter(type=readtype).order_by('sample_time')
            if len(queryset) > 0:
                stvlqs = queryset.values_list('sample_time', flat=True)
                qsvlqs = queryset.values_list('quality_sum', flat=True)
                pqsvlqs = queryset.values_list('pass_quality_sum', flat=True)
                rcvlqs = queryset.values_list('read_count', flat=True)
                prcvlqs = queryset.values_list('pass_count', flat=True)
                stcats = np.unique(stvlqs)
                unixtime = [x.timestamp() * 1000 for x in stcats]
                qs_sum = np.bincount(np.searchsorted(stcats, stvlqs), qsvlqs)
                pqs_sum = np.bincount(np.searchsorted(stcats, stvlqs), pqsvlqs)
                fqs_sum = qs_sum-pqs_sum
                rc_sum = np.bincount(np.searchsorted(stcats, stvlqs), rcvlqs)
                prc_sum = np.bincount(np.searchsorted(stcats, stvlqs), prcvlqs)
                frc_sum = rc_sum - prc_sum
                avgqs = qs_sum/rc_sum

                pavgqs = np.divide(pqs_sum, prc_sum, out=np.zeros_like(pqs_sum), where=prc_sum != 0)
                favgqs = np.divide(fqs_sum, frc_sum, out=np.zeros_like(fqs_sum), where=frc_sum != 0)
                qs = np.column_stack((unixtime,avgqs))
                pqs = np.column_stack((unixtime,pavgqs))
                fqs = np.column_stack((unixtime,favgqs))
                resultset[barcode.name]=dict()
                resultset[barcode.name][readtype.name]=dict()
                resultset[barcode.name][readtype.name]['all']=qs.tolist()
                resultset[barcode.name][readtype.name]['pass']=pqs.tolist()
                resultset[barcode.name][readtype.name]['fail']=fqs.tolist()
    return HttpResponse(json.dumps(resultset, cls=DjangoJSONEncoder), content_type="application/json")


@api_view(['GET'])
def flowcell_summary_barcode_by_minute_maxlength(request,pk):
    """
    Return a prepared summary for max read length over time grouped by minute.
    :param request:
    :param pk:
    :return:
    """
    barcodeset = BarcodeGroup.objects.filter(flowcell=pk)
    barcodedict = dict()
    resultset = dict()
    readtypes = FastqReadType.objects.all()
    for barcode in barcodeset:
        barcodedict[barcode.name] = barcode.id
        for readtype in readtypes:
            queryset = RunStatisticBarcode.objects.filter(barcode__barcodegroup=barcode.id).filter(
                type=readtype).order_by('sample_time')
            if len(queryset) > 0:
                stvlqs = queryset.values_list('sample_time', flat=True)
                maxlqs = queryset.values_list('max_length', flat=True)
                stcats = np.unique(stvlqs)
                maxarr=[]
                unixtime = [x.timestamp() * 1000 for x in stcats]
                prevtim = 0
                for maxlen, stime in zip(maxlqs, stvlqs):
                    #print (stime,prevtim)
                    if stime != prevtim:
                      maxarr.append(maxlen)
                    else:
                        #print ("got one")
                        if maxlen > maxarr[-1]:
                            maxarr[-1]=maxlen
                    prevtim = stime
                #print (len(maxarr),len(unixtime))
                qs = np.column_stack((unixtime, np.array(maxarr)))
                resultset[barcode.name] = dict()
                resultset[barcode.name][readtype.name] = dict()
                resultset[barcode.name][readtype.name]['all'] = qs.tolist()
    return HttpResponse(json.dumps(resultset, cls=DjangoJSONEncoder), content_type="application/json")


@api_view(['GET'])
def flowcell_summary_barcode_by_minute_length(request,pk):
    """
    Return a prepared set of summaries for reads quality over time grouped by minute.
    :param request:
    :param pk:
    :return:
    """
    barcodeset = BarcodeGroup.objects.filter(flowcell=pk)
    barcodedict = dict()
    resultset = dict()
    readtypes = FastqReadType.objects.all()
    for barcode in barcodeset:
        barcodedict[barcode.name] = barcode.id
        for readtype in readtypes:
            queryset = RunStatisticBarcode.objects.filter(barcode__barcodegroup=barcode.id).filter(type=readtype).order_by('sample_time')
            if len(queryset) > 0:
                stvlqs = queryset.values_list('sample_time', flat=True)
                qsvlqs = queryset.values_list('total_length', flat=True)
                pqsvlqs = queryset.values_list('pass_length', flat=True)
                rcvlqs = queryset.values_list('read_count', flat=True)
                prcvlqs = queryset.values_list('pass_count', flat=True)
                stcats = np.unique(stvlqs)
                unixtime = [x.timestamp() * 1000 for x in stcats]
                qs_sum = np.bincount(np.searchsorted(stcats, stvlqs), qsvlqs)
                pqs_sum = np.bincount(np.searchsorted(stcats, stvlqs), pqsvlqs)
                fqs_sum = qs_sum-pqs_sum
                rc_sum = np.bincount(np.searchsorted(stcats, stvlqs), rcvlqs)
                prc_sum = np.bincount(np.searchsorted(stcats, stvlqs), prcvlqs)
                frc_sum = rc_sum - prc_sum
                avgqs = qs_sum/rc_sum

                pavgqs = np.divide(pqs_sum, prc_sum, out=np.zeros_like(pqs_sum), where=prc_sum != 0)
                favgqs = np.divide(fqs_sum, frc_sum, out=np.zeros_like(fqs_sum), where=frc_sum != 0)
                qs = np.column_stack((unixtime,avgqs))
                pqs = np.column_stack((unixtime,pavgqs))
                fqs = np.column_stack((unixtime,favgqs))
                resultset[barcode.name]=dict()
                resultset[barcode.name][readtype.name]=dict()
                resultset[barcode.name][readtype.name]['all']=qs.tolist()
                resultset[barcode.name][readtype.name]['pass']=pqs.tolist()
                resultset[barcode.name][readtype.name]['fail']=fqs.tolist()
    return HttpResponse(json.dumps(resultset, cls=DjangoJSONEncoder), content_type="application/json")


@api_view(['GET'])
def flowcell_summary_barcode_by_minute_bases(request, pk):
    """
        Return a list with summaries for a particular run grouped by minute.
        """
    barcodeset = BarcodeGroup.objects.filter(flowcell=pk)
    barcodedict=dict()
    resultset=dict()
    readtypes = FastqReadType.objects.all()
    for barcode in barcodeset:
        barcodedict[barcode.name] = barcode.id
        for readtype in readtypes:
            queryset = RunStatisticBarcode.objects.filter(barcode__barcodegroup=barcode.id).filter(type=readtype).order_by('sample_time')
            if len(queryset) > 0:
                stvlqs = queryset.values_list('sample_time', flat=True)
                rcvlqs = queryset.values_list('read_count', flat=True)
                bvlqs = queryset.values_list('total_length', flat=True)
                prcvlqs = queryset.values_list('pass_count', flat=True)
                pbvlqs = queryset.values_list('pass_length', flat=True)
                stcats = np.unique(stvlqs)
                unixtime = [x.timestamp() * 1000 for x in stcats]
                rc_sum = np.bincount(np.searchsorted(stcats, stvlqs), rcvlqs)
                prc_sum = np.bincount(np.searchsorted(stcats, stvlqs), prcvlqs)
                bc_sum = np.bincount(np.searchsorted(stcats, stvlqs), bvlqs)
                pbc_sum = np.bincount(np.searchsorted(stcats, stvlqs), pbvlqs)
                fbc_sum = bc_sum - pbc_sum
                frc_sum = rc_sum - prc_sum
                rc = np.column_stack((unixtime, np.cumsum(rc_sum)))
                prc = np.column_stack((unixtime, np.cumsum(prc_sum)))
                frc = np.column_stack((unixtime, np.cumsum(frc_sum)))
                bc = np.column_stack((unixtime,np.cumsum(bc_sum)))
                pbc = np.column_stack((unixtime,np.cumsum(pbc_sum)))
                fbc = np.column_stack((unixtime,np.cumsum(fbc_sum)))
                resultset[barcode.name]=dict()
                if "reads" not in resultset[barcode.name]:
                    resultset[barcode.name]["reads"]=dict()
                resultset[barcode.name]["reads"][readtype.name]=dict()
                resultset[barcode.name]["reads"][readtype.name]['all']=rc.tolist()
                resultset[barcode.name]["reads"][readtype.name]['pass']=prc.tolist()
                resultset[barcode.name]["reads"][readtype.name]['fail']=frc.tolist()
                if "bases" not in resultset[barcode.name]:
                    resultset[barcode.name]["bases"]=dict()
                resultset[barcode.name]["bases"][readtype.name]=dict()
                resultset[barcode.name]["bases"][readtype.name]['all']=bc.tolist()
                resultset[barcode.name]["bases"][readtype.name]['pass']=pbc.tolist()
                resultset[barcode.name]["bases"][readtype.name]['fail']=fbc.tolist()

    return HttpResponse(json.dumps(resultset, cls=DjangoJSONEncoder), content_type="application/json")


@api_view(['GET'])
def flowcell_summary_barcode_by_minute_speed(request, pk):
    """
        Return a list with summaries for a particular run grouped by minute.
        """
    barcodeset = BarcodeGroup.objects.filter(flowcell=pk)
    barcodedict=dict()
    resultset=dict()
    readtypes = FastqReadType.objects.all()
    for barcode in barcodeset:
        barcodedict[barcode.name] = barcode.id
        for readtype in readtypes:
            queryset = RunStatisticBarcode.objects.filter(barcode__barcodegroup=barcode.id).filter(type=readtype).order_by('sample_time')
            if len(queryset) > 0:
                stvlqs = queryset.values_list('sample_time', flat=True)
                cpvlqs = queryset.values_list('channel_presence', flat=True)
                bvlqs = queryset.values_list('total_length', flat=True)
                #prcvlqs = queryset.values_list('pass_count', flat=True)
                pbvlqs = queryset.values_list('pass_length', flat=True)
                stcats = np.unique(stvlqs)
                channel_count = [np.count_nonzero(list(map(int,x))) for x in cpvlqs]
                #len(self.channel_presence.replace('0', ''))
                ## To increase speed, consider moving this to task and storing specific count of channels in RunStatisticBarcode
                unixtime = [x.timestamp() * 1000 for x in stcats]


                #rc_sum = np.bincount(np.searchsorted(stcats, stvlqs), rcvlqs)
                #prc_sum = np.bincount(np.searchsorted(stcats, stvlqs), prcvlqs)
                cc_sum = np.bincount(np.searchsorted(stcats,stvlqs), channel_count)
                bc_sum = np.bincount(np.searchsorted(stcats, stvlqs), bvlqs)
                pbc_sum = np.bincount(np.searchsorted(stcats, stvlqs), pbvlqs)
                fbc_sum = bc_sum - pbc_sum
                #frc_sum = rc_sum - prc_sum

                ss = np.column_stack((unixtime, bc_sum/60))
                pss = np.column_stack((unixtime, pbc_sum/60))
                fss = np.column_stack((unixtime, fbc_sum/60))
                sr = np.column_stack((unixtime, bc_sum/cc_sum/60))
                #psr = np.column_stack((unixtime, pbc_sum/cc_sum/60))
                #fsr = np.column_stack((unixtime, fbc_sum/cc_sum/60))

                resultset[barcode.name]=dict()
                if "rate" not in resultset[barcode.name]:
                    resultset[barcode.name]["rate"]=dict()
                resultset[barcode.name]["rate"][readtype.name]=dict()
                resultset[barcode.name]["rate"][readtype.name]['all']=ss.tolist()
                resultset[barcode.name]["rate"][readtype.name]['pass']=pss.tolist()
                resultset[barcode.name]["rate"][readtype.name]['fail']=fss.tolist()
                if "speed" not in resultset[barcode.name]:
                    resultset[barcode.name]["speed"]=dict()
                resultset[barcode.name]["speed"][readtype.name]=dict()
                resultset[barcode.name]["speed"][readtype.name]['all']=sr.tolist()
                #resultset[barcode.name]["speed"][readtype.name]['pass']=psr.tolist()
                #resultset[barcode.name]["speed"][readtype.name]['fail']=fsr.tolist()

    return HttpResponse(json.dumps(resultset, cls=DjangoJSONEncoder), content_type="application/json")


@api_view(['GET'])
def flowcell_histogram_summary(request, pk):
    """
    Return a list with histogram summaries for a particular run.
    """
    queryset = FlowCellRun.objects.filter(flowcell_id=pk)
    runset = list()
    for run in queryset:
        # print (run.run_id)
        runset.append(run.run_id)
    queryset = HistogramSummary.objects\
        .filter(run_id__owner=request.user)\
        .filter(run_id__in=runset)\
        .order_by('read_type', 'bin_width')

    serializer = RunHistogramSummarySerializer(queryset, many=True, context={'request': request})

    return Response(serializer.data)

@api_view(['GET'])
def flowcell_channel_summary(request, pk):
    """
    Return a list with channel info for a particular run.
    """
    queryset = FlowCellRun.objects.filter(flowcell_id=pk)
    runset = list()
    for run in queryset:
        # print (run.run_id)
        runset.append(run.run_id)
    queryset = ChannelSummary.objects\
        .filter(run_id__owner=request.user)\
        .filter(run_id__in=runset)

    serializer = ChannelSummarySerializer(queryset, many=True, context={'request': request})

    return Response(serializer.data)



@api_view(['GET'])
def flowcell_run_status_list(request,pk):
    """
    TODO describe function

    """
    queryset = FlowCellRun.objects.filter(flowcell_id=pk)
    runset = list()
    for run in queryset:
        # print (run.run_id)
        runset.append(run.run_id)

    try:
        minIONrunstat = MinIONRunStatus.objects.filter(run_id__in=runset)

    except MinIONRunStatus.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        serializer = MinIONRunStatusSerializer(minIONrunstat, many=True, context={'request': request})
        return Response(serializer.data)



@api_view(['GET'])
def flowcell_channel_summary_readcount(request, pk):
    """
    Return a list with channel info for a particular run.
    """
    queryset = FlowCellRun.objects.filter(flowcell_id=pk)
    runset = list()
    flowcell_type = 0
    for run in queryset:
        runset.append(run.run_id)
        if run.run.flowcell_type() > flowcell_type:
            flowcell_type = run.run.flowcell_type()
    queryset = ChannelSummary.objects\
        .filter(run_id__owner=request.user)\
        .filter(run_id__in=runset)
    chan_events=dict()
    query_value_list=queryset.values_list('channel_number',flat=True)
    for query in queryset:
        if query.channel_number not in chan_events.keys():
            chan_events[query.channel_number]=dict()
            chan_events[query.channel_number]["read_count"]=query.read_count
        else:
            chan_events[query.channel_number]["read_count"] += query.read_count
    datatoreturn=list()
    for i in range(flowcell_type):
        if (i+1) in chan_events:
            coordinate = get_coords((i+1),flowcell_type)
            datatoreturn.append((coordinate[0],coordinate[1],chan_events[i + 1]['read_count']))
        else:
            coordinate = get_coords((i+1),flowcell_type)
            datatoreturn.append((coordinate[0], coordinate[1], 0))
    return HttpResponse(json.dumps(datatoreturn),content_type="application/json")


@api_view(['GET'])
def flowcell_channel_summary_readkb(request, pk):
    """
    Return a list with channel info for a particular run.
    """
    queryset = FlowCellRun.objects.filter(flowcell_id=pk)
    runset = list()
    flowcell_type = 0
    for run in queryset:
        runset.append(run.run_id)
        if run.run.flowcell_type() > flowcell_type:
            flowcell_type = run.run.flowcell_type()
    queryset = ChannelSummary.objects\
        .filter(run_id__owner=request.user)\
        .filter(run_id__in=runset)
    chan_events=dict()
    query_value_list=queryset.values_list('channel_number',flat=True)
    for query in queryset:
        if query.channel_number not in chan_events.keys():
            chan_events[query.channel_number]=dict()
            chan_events[query.channel_number]["read_length"]=round(query.read_length/1000)
        else:
            chan_events[query.channel_number]["read_length"] += round(query.read_length / 1000)
    datatoreturn=list()
    for i in range(flowcell_type):
        if (i+1) in chan_events:
            coordinate = get_coords((i+1),flowcell_type)
            datatoreturn.append((coordinate[0],coordinate[1],chan_events[i + 1]['read_length']))
        else:
            coordinate = get_coords((i+1),flowcell_type)
            datatoreturn.append((coordinate[0], coordinate[1], 0))
    return HttpResponse(json.dumps(datatoreturn),content_type="application/json")

def get_coords(channel,flowcellsize):
    #if flowcellsize==3000:
    #    return (channel//30, channel%30)
    if flowcellsize==3000:
        #find which block of 12 we are in:
        block=(channel-1)//250
        remainder=(channel-1)%250
        row=remainder//10
        column=remainder%10+block*10
        return(column,row)
    if flowcellsize==128:
        return (channel//8, channel%8)
    else:
        return minION_flowcell_layout(channel)




def minION_flowcell_layout(channel):
    chanlookup = {1: (31, 0), 2: (31, 1), 3: (31, 2), 4: (31, 3), 5: (31, 4), 6: (31, 5), 7: (31, 6), 8: (31, 7),
                  9: (30, 0), 10: (30, 1), 11: (30, 2), 12: (30, 3), 13: (30, 4), 14: (30, 5), 15: (30, 6), 16: (30, 7),
                  17: (29, 0), 18: (29, 1), 19: (29, 2), 20: (29, 3), 21: (29, 4), 22: (29, 5), 23: (29, 6),
                  24: (29, 7), 25: (28, 0), 26: (28, 1), 27: (28, 2), 28: (28, 3), 29: (28, 4), 30: (28, 5),
                  31: (28, 6), 32: (28, 7), 33: (31, 15), 34: (31, 14), 35: (31, 13), 36: (31, 12), 37: (31, 11),
                  38: (31, 10), 39: (31, 9), 40: (31, 8), 41: (30, 15), 42: (30, 14), 43: (30, 13), 44: (30, 12),
                  45: (30, 11), 46: (30, 10), 47: (30, 9), 48: (30, 8), 49: (29, 15), 50: (29, 14), 51: (29, 13),
                  52: (29, 12), 53: (29, 11), 54: (29, 10), 55: (29, 9), 56: (29, 8), 57: (28, 15), 58: (28, 14),
                  59: (28, 13), 60: (28, 12), 61: (28, 11), 62: (28, 10), 63: (28, 9), 64: (28, 8), 65: (3, 0),
                  66: (3, 1), 67: (3, 2), 68: (3, 3), 69: (3, 4), 70: (3, 5), 71: (3, 6), 72: (3, 7), 73: (2, 0),
                  74: (2, 1), 75: (2, 2), 76: (2, 3), 77: (2, 4), 78: (2, 5), 79: (2, 6), 80: (2, 7), 81: (1, 0),
                  82: (1, 1), 83: (1, 2), 84: (1, 3), 85: (1, 4), 86: (1, 5), 87: (1, 6), 88: (1, 7), 89: (0, 0),
                  90: (0, 1), 91: (0, 2), 92: (0, 3), 93: (0, 4), 94: (0, 5), 95: (0, 6), 96: (0, 7), 97: (3, 15),
                  98: (3, 14), 99: (3, 13), 100: (3, 12), 101: (3, 11), 102: (3, 10), 103: (3, 9), 104: (3, 8),
                  105: (2, 15), 106: (2, 14), 107: (2, 13), 108: (2, 12), 109: (2, 11), 110: (2, 10), 111: (2, 9),
                  112: (2, 8), 113: (1, 15), 114: (1, 14), 115: (1, 13), 116: (1, 12), 117: (1, 11), 118: (1, 10),
                  119: (1, 9), 120: (1, 8), 121: (0, 15), 122: (0, 14), 123: (0, 13), 124: (0, 12), 125: (0, 11),
                  126: (0, 10), 127: (0, 9), 128: (0, 8), 129: (7, 0), 130: (7, 1), 131: (7, 2), 132: (7, 3),
                  133: (7, 4), 134: (7, 5), 135: (7, 6), 136: (7, 7), 137: (6, 0), 138: (6, 1), 139: (6, 2),
                  140: (6, 3), 141: (6, 4), 142: (6, 5), 143: (6, 6), 144: (6, 7), 145: (5, 0), 146: (5, 1),
                  147: (5, 2), 148: (5, 3), 149: (5, 4), 150: (5, 5), 151: (5, 6), 152: (5, 7), 153: (4, 0),
                  154: (4, 1), 155: (4, 2), 156: (4, 3), 157: (4, 4), 158: (4, 5), 159: (4, 6), 160: (4, 7),
                  161: (7, 15), 162: (7, 14), 163: (7, 13), 164: (7, 12), 165: (7, 11), 166: (7, 10), 167: (7, 9),
                  168: (7, 8), 169: (6, 15), 170: (6, 14), 171: (6, 13), 172: (6, 12), 173: (6, 11), 174: (6, 10),
                  175: (6, 9), 176: (6, 8), 177: (5, 15), 178: (5, 14), 179: (5, 13), 180: (5, 12), 181: (5, 11),
                  182: (5, 10), 183: (5, 9), 184: (5, 8), 185: (4, 15), 186: (4, 14), 187: (4, 13), 188: (4, 12),
                  189: (4, 11), 190: (4, 10), 191: (4, 9), 192: (4, 8), 193: (11, 0), 194: (11, 1), 195: (11, 2),
                  196: (11, 3), 197: (11, 4), 198: (11, 5), 199: (11, 6), 200: (11, 7), 201: (10, 0), 202: (10, 1),
                  203: (10, 2), 204: (10, 3), 205: (10, 4), 206: (10, 5), 207: (10, 6), 208: (10, 7), 209: (9, 0),
                  210: (9, 1), 211: (9, 2), 212: (9, 3), 213: (9, 4), 214: (9, 5), 215: (9, 6), 216: (9, 7),
                  217: (8, 0), 218: (8, 1), 219: (8, 2), 220: (8, 3), 221: (8, 4), 222: (8, 5), 223: (8, 6),
                  224: (8, 7), 225: (11, 15), 226: (11, 14), 227: (11, 13), 228: (11, 12), 229: (11, 11), 230: (11, 10),
                  231: (11, 9), 232: (11, 8), 233: (10, 15), 234: (10, 14), 235: (10, 13), 236: (10, 12), 237: (10, 11),
                  238: (10, 10), 239: (10, 9), 240: (10, 8), 241: (9, 15), 242: (9, 14), 243: (9, 13), 244: (9, 12),
                  245: (9, 11), 246: (9, 10), 247: (9, 9), 248: (9, 8), 249: (8, 15), 250: (8, 14), 251: (8, 13),
                  252: (8, 12), 253: (8, 11), 254: (8, 10), 255: (8, 9), 256: (8, 8), 257: (15, 0), 258: (15, 1),
                  259: (15, 2), 260: (15, 3), 261: (15, 4), 262: (15, 5), 263: (15, 6), 264: (15, 7), 265: (14, 0),
                  266: (14, 1), 267: (14, 2), 268: (14, 3), 269: (14, 4), 270: (14, 5), 271: (14, 6), 272: (14, 7),
                  273: (13, 0), 274: (13, 1), 275: (13, 2), 276: (13, 3), 277: (13, 4), 278: (13, 5), 279: (13, 6),
                  280: (13, 7), 281: (12, 0), 282: (12, 1), 283: (12, 2), 284: (12, 3), 285: (12, 4), 286: (12, 5),
                  287: (12, 6), 288: (12, 7), 289: (15, 15), 290: (15, 14), 291: (15, 13), 292: (15, 12), 293: (15, 11),
                  294: (15, 10), 295: (15, 9), 296: (15, 8), 297: (14, 15), 298: (14, 14), 299: (14, 13), 300: (14, 12),
                  301: (14, 11), 302: (14, 10), 303: (14, 9), 304: (14, 8), 305: (13, 15), 306: (13, 14), 307: (13, 13),
                  308: (13, 12), 309: (13, 11), 310: (13, 10), 311: (13, 9), 312: (13, 8), 313: (12, 15), 314: (12, 14),
                  315: (12, 13), 316: (12, 12), 317: (12, 11), 318: (12, 10), 319: (12, 9), 320: (12, 8), 321: (19, 0),
                  322: (19, 1), 323: (19, 2), 324: (19, 3), 325: (19, 4), 326: (19, 5), 327: (19, 6), 328: (19, 7),
                  329: (18, 0), 330: (18, 1), 331: (18, 2), 332: (18, 3), 333: (18, 4), 334: (18, 5), 335: (18, 6),
                  336: (18, 7), 337: (17, 0), 338: (17, 1), 339: (17, 2), 340: (17, 3), 341: (17, 4), 342: (17, 5),
                  343: (17, 6), 344: (17, 7), 345: (16, 0), 346: (16, 1), 347: (16, 2), 348: (16, 3), 349: (16, 4),
                  350: (16, 5), 351: (16, 6), 352: (16, 7), 353: (19, 15), 354: (19, 14), 355: (19, 13), 356: (19, 12),
                  357: (19, 11), 358: (19, 10), 359: (19, 9), 360: (19, 8), 361: (18, 15), 362: (18, 14), 363: (18, 13),
                  364: (18, 12), 365: (18, 11), 366: (18, 10), 367: (18, 9), 368: (18, 8), 369: (17, 15), 370: (17, 14),
                  371: (17, 13), 372: (17, 12), 373: (17, 11), 374: (17, 10), 375: (17, 9), 376: (17, 8), 377: (16, 15),
                  378: (16, 14), 379: (16, 13), 380: (16, 12), 381: (16, 11), 382: (16, 10), 383: (16, 9), 384: (16, 8),
                  385: (23, 0), 386: (23, 1), 387: (23, 2), 388: (23, 3), 389: (23, 4), 390: (23, 5), 391: (23, 6),
                  392: (23, 7), 393: (22, 0), 394: (22, 1), 395: (22, 2), 396: (22, 3), 397: (22, 4), 398: (22, 5),
                  399: (22, 6), 400: (22, 7), 401: (21, 0), 402: (21, 1), 403: (21, 2), 404: (21, 3), 405: (21, 4),
                  406: (21, 5), 407: (21, 6), 408: (21, 7), 409: (20, 0), 410: (20, 1), 411: (20, 2), 412: (20, 3),
                  413: (20, 4), 414: (20, 5), 415: (20, 6), 416: (20, 7), 417: (23, 15), 418: (23, 14), 419: (23, 13),
                  420: (23, 12), 421: (23, 11), 422: (23, 10), 423: (23, 9), 424: (23, 8), 425: (22, 15), 426: (22, 14),
                  427: (22, 13), 428: (22, 12), 429: (22, 11), 430: (22, 10), 431: (22, 9), 432: (22, 8), 433: (21, 15),
                  434: (21, 14), 435: (21, 13), 436: (21, 12), 437: (21, 11), 438: (21, 10), 439: (21, 9), 440: (21, 8),
                  441: (20, 15), 442: (20, 14), 443: (20, 13), 444: (20, 12), 445: (20, 11), 446: (20, 10),
                  447: (20, 9), 448: (20, 8), 449: (27, 0), 450: (27, 1), 451: (27, 2), 452: (27, 3), 453: (27, 4),
                  454: (27, 5), 455: (27, 6), 456: (27, 7), 457: (26, 0), 458: (26, 1), 459: (26, 2), 460: (26, 3),
                  461: (26, 4), 462: (26, 5), 463: (26, 6), 464: (26, 7), 465: (25, 0), 466: (25, 1), 467: (25, 2),
                  468: (25, 3), 469: (25, 4), 470: (25, 5), 471: (25, 6), 472: (25, 7), 473: (24, 0), 474: (24, 1),
                  475: (24, 2), 476: (24, 3), 477: (24, 4), 478: (24, 5), 479: (24, 6), 480: (24, 7), 481: (27, 15),
                  482: (27, 14), 483: (27, 13), 484: (27, 12), 485: (27, 11), 486: (27, 10), 487: (27, 9), 488: (27, 8),
                  489: (26, 15), 490: (26, 14), 491: (26, 13), 492: (26, 12), 493: (26, 11), 494: (26, 10),
                  495: (26, 9), 496: (26, 8), 497: (25, 15), 498: (25, 14), 499: (25, 13), 500: (25, 12), 501: (25, 11),
                  502: (25, 10), 503: (25, 9), 504: (25, 8), 505: (24, 15), 506: (24, 14), 507: (24, 13), 508: (24, 12),
                  509: (24, 11), 510: (24, 10), 511: (24, 9), 512: (24, 8)}
    return chanlookup[channel]

@api_view(['GET'])
def flowcell_run_status_list(request, pk):
    if request.method == 'GET':
        queryset = FlowCellRun.objects.filter(flowcell_id=pk)
        runset = list()
        for run in queryset:
            # print (run.run_id)
            runset.append(run.run_id)
        queryset = MinIONRunStatus.objects.filter(run_id__in=runset)
        #print (queryset)
        serializer = MinIONRunStatusSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)


@api_view(['GET'])
def flowcell_run_stats_list(request,pk):
    """
    TODO describe function
    """
    queryset = FlowCellRun.objects.filter(flowcell_id=pk)
    runset = list()
    for run in queryset:
        # print (run.run_id)
        runset.append(run.run_id)

    try:
        crazyminIONrunstats = MinIONRunStats.objects.filter(run_id__in=runset)
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


@api_view(['GET'])
def flowcell_run_stats_latest(request,pk,checkid):
    """
    TODO describe function
    """
    queryset = FlowCellRun.objects.filter(flowcell_id=pk)
    runset = list()
    for run in queryset:
        # print (run.run_id)
        runset.append(run.run_id)
    try:
        crazyminIONrunstats = MinIONRunStats.objects.filter(run_id__in=runset, id__gt=checkid)[:1000]
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
        flowcellitem=FlowCell.objects.get(id=pk)
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
def flowcell_tasks_detail_all(request,pk):
    queryset = FlowCellRun.objects.filter(flowcell_id=pk)
    runset = list()
    for run in queryset:
        # print (run.run_id)
        runset.append(run.run_id)

    queryset=JobType.objects.filter(private=False).filter(type_name__type_name__in=['flowcell',])
    minionruns=MinIONRun.objects.filter(run_id__in=runset)

    result = []
    #print (queryset)
    for jobtype in queryset:
        #print (jobtype.type_name.all())
        obj = {}
        obj.update({
            'name': jobtype.name,
            'description' : jobtype.description,
            'long_description' : jobtype.long_description,
            'read_count': jobtype.readcount,
            'reference': jobtype.reference,
            'transcriptome': jobtype.transcriptome
        })

        jobmasterlist = JobMaster.objects.filter(Q(run_id__in=runset) | Q(flowcell_id=pk)).filter(job_type=jobtype)

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
        }
    }
    tabs = list()
    tabs_send = list()
    queryset = FlowCellRun.objects.filter(flowcell_id=pk)
    runset = list()

    for run in queryset:
        runset.append(run.run_id)

    if MinIONRunStatus.objects.filter(run_id__in=runset):
        tabs.append(flowcell_tabs_dict['LiveEvent'])

    for master in JobMaster.objects.filter(Q(run_id__in=runset) | Q(flowcell_id=pk)).filter(last_read__gt=0).values_list('job_type__name', flat=True):
        if master in flowcell_tabs_dict.keys():
            tabs.append(flowcell_tabs_dict[master])
        else:
            print("Flowcell '" + pk + "' has JobType '" + master + "' but there is no corresponding tab defined in reads/views.py")

    for tab in tabs:
        if tab not in tabs_send:
            tabs_send.append(tab)

    return Response(tabs_send)
