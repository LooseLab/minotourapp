from django.shortcuts import render
from django.db.models import Count, Min, Sum, Avg

# Create your views here.

from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from alignment.models import PafRoughCov
from alignment.models import PafSummaryCov
from alignment.serializers import PafRoughCovSerializer
from alignment.serializers import PafRoughCovChromSerializer
from alignment.serializers import PafRoughCovChromSerializerCount
from alignment.serializers import PafSummaryCovSerializer

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
def paf_alignment_summary(request, pk):#,bc,ch):
    """

    :param request:
    :param pk:
    :return:
    """
    if request.method == 'GET':
        queryset = PafSummaryCov.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk)
            #.filter(run__id=pk,barcode=bc,chromosome=ch)

        serializer = PafSummaryCovSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)

@api_view(['GET'])
def paf_test(request, pk,bc,ch,ty,po,ln):
    """

    :param request:
    :param pk:
    :return:
    """
    if request.method == 'GET':
        queryset = PafRoughCov.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk) \
            .filter(barcode__id=bc) \
            .filter(chromosome__id=ch) \
            .filter(read_type__id=ty) \
            .filter(p__lte=po).values('reference').annotate(sumcount=Sum('i'))
        # .filter(run__id=pk,barcode=bc,chromosome=ch)

        #serializer = PafRoughCovSerializer(queryset, many=True, context={'request': request})

        #queryset = PafRoughCov.objects.filter(run__owner=request.user).filter(run__id=pk).filter(p__lte=1000000).values('reference').annotate(sumcount = Sum('i')).order_by('sumcount')

        serializer1 = PafRoughCovChromSerializerCount(queryset, many=True, context={'request': request})

        queryset2 = PafRoughCov.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk) \
            .filter(barcode__id=bc) \
            .filter(chromosome__id=ch) \
            .filter(read_type__id=ty) \
            .filter(p__gt=po)\
            .filter(p__lt=(po+ln))\
            .order_by('p')

        serializer2 = PafRoughCovChromSerializer(queryset2, many=True, context={'request': request})

        Serializer_list = [serializer1.data, serializer2.data]

        content = {
            'status': 1,
            'responseCode': status.HTTP_200_OK,
            'data': Serializer_list,

        }

        return Response(content)


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