from django.shortcuts import render

# Create your views here.

from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from alignment.models import PafRoughCov
from alignment.serializers import PafRoughCovSerializer
from alignment.serializers import PafRoughCovChromSerializer

from reference.models import ReferenceInfo
from reference.models import ReferenceLine

@api_view(['GET'])
def paf_alignment_list(request, pk):#,bc,ch):
    """

    :param request:
    :param pk:
    :return:
    """
    if request.method == 'GET':
        queryset = PafRoughCov.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk)
            #.filter(run__id=pk,barcode=bc,chromosome=ch)

        serializer = PafRoughCovSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)


@api_view(['GET'])
def paf_alignment_ch_list(request, pk,ch):#,bc,ch):
    """

    :param request:
    :param pk:
    :return:
    """
    if request.method == 'GET':
        queryset = PafRoughCov.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk).filter(chromosome__id=ch) \
            .order_by('p')
            #.filter(run__id=pk,barcode=bc,chromosome=ch)

        serializer = PafRoughCovChromSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)