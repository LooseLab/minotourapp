from django.shortcuts import render
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from reference.models import ReferenceInfo, ReferenceLine
from reference.serializers import (ReferenceInfoSerializer,
                                   ReferenceLineSerializer)

# Create your views here.



@api_view(['GET'])
def reference_list(request, pk):#,bc,ch):
    """

    :param request:
    :param pk:
    :return:
    """
    if request.method == 'GET':
        queryset = ReferenceInfo.objects \
            .filter(id=pk)

        serializer = ReferenceInfoSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)


@api_view(['GET'])
def reference_line(request, pk):
    """

    :param request:
    :param pk:
    :return:
    """
    if request.method == 'GET':
        queryset = ReferenceLine.objects.filter(pk=pk)
        serializer = ReferenceLineSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)

@api_view(['GET'])
def reference_all(request):
    if request.method == 'GET':
        queryset = ReferenceInfo.objects
        serializer = ReferenceInfoSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)
