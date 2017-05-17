from django.shortcuts import render
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import authentication
from rest_framework import permissions
from rest_framework import status
from rest_framework import viewsets

from reads.models import MinionRun
from reads.serializers import MinionRunSerializer


class RunsList(APIView):
    """
    Returns a list of runs (GET) or create a new run (POST).
    """

    authentication_classes = (authentication.TokenAuthentication,)
    permission_classes = (permissions.IsAuthenticated,)

    def get(self, request):
        """
        Return a list of all run of the user associated with the token.
        """

        queryset = MinionRun.objects.filter(owner=request.user)

        serializer = MinionRunSerializer(queryset, many=True)

        return Response(serializer.data)

    def post(self, request):
        """
        Create a new run which owner is the user associated with the token.
        """

        serializer = MinionRunSerializer(data=request.data)

        if serializer.is_valid():
            serializer.save(owner=request.user)
            return Response(serializer.data, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
