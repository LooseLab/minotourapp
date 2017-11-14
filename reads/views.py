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

from reads.models import Barcode, ChannelSummary
from reads.models import FastqRead
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
from reads.serializers import BarcodeSerializer
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
#def read_list(request, pk,page):
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

        print('Serializer errors: {} '.format(serializer.errors))
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
        serializer = BarcodeSerializer(data=request.data, context={'request': request})

        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


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
    queryset=JobType.objects.filter(private=False)
    minionrun=MinIONRun.objects.get(pk=pk)

    result = []

    for jobtype in queryset:
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
                'temp_file': jobmasterlist[0].temp_file,
                'complete': jobmasterlist[0].complete,
                'running': jobmasterlist[0].running
            })

            obj.update({
                'job_details': obj2
            })

        result.append(obj)

    return HttpResponse(json.dumps(result), content_type="application/json")


@api_view(['GET','POST'])
def flowcell_list(request):
    if request.method == 'GET':
        queryset = FlowCell.objects.filter(owner=request.user)
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


    #alldata = MinIONRun.objects.filter(id__in=runset)

    #serializer = MinIONRunSerializer(alldata,many=True, context={'request': request})
    #return Response(serializer.data)


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
def flowcell_run_status_list(request, pk):
    if request.method == 'GET':
        queryset = FlowCellRun.objects.filter(flowcell_id=pk)
        runset = list()
        for run in queryset:
            # print (run.run_id)
            runset.append(run.run_id)
        queryset = MinIONRunStatus.objects.filter(run_id__in=runset)
        print (queryset)
        serializer = MinIONRunStatusSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)