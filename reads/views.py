import json
from dateutil import parser


from django.http import HttpResponse
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from reads.models import MinionRun, FastqRead, FastqReadType, RunStatistic
from reads.serializers import MinionRunSerializer, FastqReadSerializer, FastqReadTypeSerializer


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
        queryset = MinionRun.objects.filter(owner=request.user)
        serializer = MinionRunSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'POST':
        serializer = MinionRunSerializer(data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save(owner=request.user)
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(['GET', 'PUT', 'DELETE'])
def run_detail(request, pk):
    """
    Retrieve, update or delete a run instance.
    """
    try:
        run = MinionRun.objects.get(pk=pk)

        if run.owner != request.user:
            return Response(status=status.HTTP_401_UNAUTHORIZED)

    except MinionRun.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        serializer = MinionRunSerializer(run, context={'request': request})
        return Response(serializer.data)

    elif request.method == 'PUT':
        serializer = MinionRunSerializer(run, data=request.data, context={'request': request})
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    elif request.method == 'DELETE':
        run.delete()
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
    if request.method == 'GET':
        queryset = FastqRead.objects.filter(run_id=pk)
        result=set()
        for key in queryset:
            result.add(key.read_id)
        return HttpResponse(json.dumps(list(result)), content_type="application/json")

@api_view(['GET'])
def cumulative_read_count(request,pk):
    if request.method == 'GET':
        queryset = RunStatistic.objects.filter(run_id__owner=request.user).filter(run_id=pk).order_by('type','sample_time',)
        result=dict()
        for key in queryset:
            if str(key.type) not in result:
                result[str(key.type)]=dict()
            result[str(key.type)][str(key.sample_time)]=key.total_length
        data_to_return=dict()
        counter = 0
        for readtype in result:
            cumuyield=0
            data_to_return[counter]=dict()
            data_to_return[counter]['name']=readtype
            data_to_return[counter]['data']=list()
            for data in result[readtype]:
                cumuyield = cumuyield + result[readtype][data]
                data_to_return[counter]['data'].append((UTC_time_to_epoch(data), cumuyield))
        return HttpResponse(json.dumps(data_to_return), content_type="application/json")

def UTC_time_to_epoch(timestamp):
    dt = parser.parse(timestamp)
    return dt.timestamp()*1000