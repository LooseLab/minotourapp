from django.shortcuts import render
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response
# Create your views here.

from minikraken.models import MiniKraken
from minikraken.models import ParsedKraken

from minikraken.serializers import MiniKrakenSerializer
from minikraken.serializers import ParsedKrakenSerializer


@api_view(['GET'])
def kraken_raw_list(request, pk):#,bc,ch):
    """

    :param request:
    :param pk:
    :return:
    """
    if request.method == 'GET':
        queryset = MiniKraken.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk)
            #.filter(run__id=pk,barcode=bc,chromosome=ch)

        serializer = MiniKrakenSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)


@api_view(['GET'])
def kraken_parse_list(request, pk):#,bc,ch):
    """

    :param request:
    :param pk:
    :return:
    """
    if request.method == 'GET':
        queryset = ParsedKraken.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk)
            #.filter(run__id=pk,barcode=bc,chromosome=ch)

        serializer = ParsedKrakenSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)



@api_view(['GET'])
def flowcellkraken_raw_list(request, pk):#,bc,ch):
    """

    :param request:
    :param pk:
    :return:
    """


    if request.method == 'GET':
        queryset = MiniKraken.objects \
            .filter(run__owner=request.user) \
            .filter(flowcell__id=pk)
            #.filter(run__id=pk,barcode=bc,chromosome=ch)

        serializer = MiniKrakenSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)


@api_view(['GET'])
def flowcellkraken_parse_list(request, pk):#,bc,ch):
    """

    :param request:
    :param pk:
    :return:
    """
    if request.method == 'GET':
        queryset = ParsedKraken.objects \
            .filter(flowcell__id=pk) \
            .filter(percentage__gte=0.03)
            #.filter(run__id=pk,barcode=bc,chromosome=ch)

        serializer = ParsedKrakenSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)

