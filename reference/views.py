from django.db.models import Q
from rest_framework.decorators import api_view
from rest_framework.response import Response

from reference.models import ReferenceInfo, ReferenceLine
from reference.serializers import (ReferenceInfoSerializer,
                                   ReferenceLineSerializer)


@api_view(['GET'])
def reference_list(request):
    """
    Get a list of all references available to the User making the request.
    Parameters
    ----------
    request: rest_framework.request.Request
        Query params - names. Whether we only want the names and IDs of the references.
    Returns
    -------
    rest_framework.response.Response
    """

    if request.GET.get("names", False):
        data = ReferenceInfo.objects.filter(Q(private=False) | Q(owner=request.user)).values_list("name", "id")
    else:
        queryset = ReferenceInfo.objects.filter(Q(private=False) | Q(owner=request.user))

        serializer = ReferenceInfoSerializer(queryset, many=True, context={'request': request})

        data = serializer.data

    return Response(data)


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

