import json

from django.http import HttpResponse
from django.shortcuts import render
from django.db.models import Count, Min, Sum, Avg

from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from assembly.models import GfaSummary
from assembly.serializers import GfaSummarySerializer
from assembly.serializers import FlowcellGfaSummarySerializer


@api_view(['GET'])
def gfa_output_list(request, run_id):

    queryset = GfaSummary.objects \
        .filter(run__owner=request.user) \
        .filter(run__id=run_id) \
        .order_by('nreads')

    # ncontigs_list = []
    # position_list = []
    # coverage_list = []
    # current_coverage_sum = 0
    #
    # for key, item in enumerate(queryset):
    #
    #     ncontigs_list.append([item.nreads, item.ncontigs])

    # return HttpResponse(json.dumps(ncontigs_list), content_type="application/json")

    serializer = GfaSummarySerializer(queryset, many=True, context={'request': request})

    return Response(serializer.data)


@api_view(['GET'])
def flowcell_gfa_output_list(request, run_id):

    queryset = GfaSummary.objects \
        .filter(flowcell__id=run_id) \
        .order_by('nreads')

    serializer = FlowcellGfaSummarySerializer(queryset, many=True, context={'request': request})
    print (serializer.data)
    print (run_id)
    print (request.user)
    print (queryset)

    return Response(serializer.data)
