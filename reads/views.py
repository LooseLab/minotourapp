from rest_framework import generics
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response

from reads.models import MinionRun, FastqRead
from reads.permissions import IsOwner
from reads.serializers import MinionRunSerializer, FastqReadSerializer


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


class ReadListView(generics.ListCreateAPIView):
    queryset = FastqRead.objects.all()
    serializer_class = FastqRead
    permission_classes = (IsAuthenticated,)

    def list(self, request, *args, **kwargs):
        queryset = FastqRead.objects.filter(run_id__id=self.kwargs['pk'])
        serializer = FastqReadSerializer(queryset, many=True)
        return Response(serializer.data)
