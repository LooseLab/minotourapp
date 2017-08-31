import json
from datetime import datetime, timedelta
from django.utils import timezone
from dateutil import parser
from django.db.models import Q
from django.http import HttpResponse
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from reads.models import FastqRead, RunSummaryBarCode, RunStatisticBarcode
from reads.models import FastqReadType
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
from reads.models import RunSummary
from reads.serializers import FastqReadSerializer, RunSummaryBarcodeSerializer, RunStatisticBarcodeSerializer, \
    RunStatisticSerializer
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
from reads.serializers import RunSummarySerializer



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
        queryset = MinIONRun.objects.filter(owner=request.user)
        serializer = MinIONRunSerializer(queryset, many=True, context={'request': request})
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
        queryset = MinIONRun.objects.filter(owner=request.user).filter(active=True).distinct()
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
        queryset = MinION.objects.filter(owner=request.user).exclude(events__event__name='unplugged')
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


@api_view(['GET', 'PUT', 'DELETE'])
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
        queryset = FastqRead.objects.filter(run_id=pk)
        serializer = FastqReadSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'POST':
        serializer = FastqReadSerializer(data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
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
    Retrieve, update or delete a run instance.
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


@api_view(['GET'])
def readname_list(request, pk):
    """
    TODO describe function
    """
    if request.method == 'GET':
        queryset = FastqRead.objects.filter(run_id=pk)
        result = set()
        for key in queryset:
            result.add(key.read_id)
        return HttpResponse(json.dumps(list(result)), content_type="application/json")


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
    queryset = RunSummaryBarCode.objects\
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
