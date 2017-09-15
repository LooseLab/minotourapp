from django.shortcuts import render

# Create your views here.

from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from alignment.models import PafRoughCov
from alignment.serializers import PafRoughCovSerializer

@api_view(['GET'])
def paf_alignment_list(request, pk):
    """

    :param request:
    :param pk:
    :return:
    """
    if request.method == 'GET':
        queryset = PafRoughCov.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk)

        serializer = PafRoughCovSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)