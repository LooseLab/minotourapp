from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from readuntil.serializers import RejectedFastqReadSerializer


@api_view(['POST'])
def rejectedfastqread_list(request):

    serializer = RejectedFastqReadSerializer(
        data=request.data,
        context={'request': request}
    )

    if serializer.is_valid():

        serializer.save(owner=request.user)
        return Response(serializer.data, status=status.HTTP_201_CREATED)

    return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


