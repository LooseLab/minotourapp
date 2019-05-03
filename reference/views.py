from django.db.models import Q
from rest_framework.decorators import api_view
from rest_framework.response import Response

from reference.models import ReferenceInfo, ReferenceLine
from reference.serializers import (ReferenceInfoSerializer,
                                   ReferenceLineSerializer)


@api_view(['GET'])
def reference_list(request):

    queryset = ReferenceInfo.objects.filter(Q(private=False) | Q(owner=request.user.id))

    serializer = ReferenceInfoSerializer(queryset, many=True, context={'request': request})

    return Response(serializer.data)


@api_view(['GET'])
def reference_detail(request, pk):

    queryset = ReferenceInfo.objects.filter(Q(private=False) | Q(owner=request.user.id), pk=pk)[0]

    serializer = ReferenceInfoSerializer(queryset, context={'request': request})

    return Response(serializer.data)


@api_view(['GET'])
def reference_line(request, pk):

    if request.method == 'GET':

        queryset = ReferenceLine.objects.filter(pk=pk)

        serializer = ReferenceLineSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)

