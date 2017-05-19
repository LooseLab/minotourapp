from rest_framework import generics
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response

from reads.models import MinionRun
from reads.permissions import IsOwner
from reads.serializers import MinionRunSerializer


class RunListView(generics.ListCreateAPIView):
    queryset = MinionRun.objects.all()
    serializer_class = MinionRunSerializer
    permission_classes = (IsAuthenticated,)

    def list(self, request, *args, **kwargs):
        queryset = MinionRun.objects.filter(owner=request.user)
        serializer = MinionRunSerializer(queryset, many=True)
        return Response(serializer.data)

    def perform_create(self, serializer):
        serializer.save(owner=self.request.user)


class RunDetailView(generics.RetrieveUpdateDestroyAPIView):
    queryset = MinionRun.objects.all()
    serializer_class = MinionRunSerializer
    permission_classes = (IsOwner,)
