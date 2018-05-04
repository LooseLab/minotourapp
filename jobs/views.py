from rest_framework.decorators import api_view
from rest_framework.response import Response

from jobs.models import JobType
from reads.serializers import JobTypeSerializer


@api_view(['GET'])
def job_type_list(request):

    if request.method == 'GET':

        job_type_list = JobType.objects.filter(private=False)

        serializer = JobTypeSerializer(
            job_type_list,
            many=True,
            context={
                'request': request
            }
        )

        return Response(serializer.data)


@api_view(['GET'])
def job_type_detail(request, id):

    if request.method == 'GET':

        job_type_list = JobType.objects.filter(private=False)

        serializer = JobTypeSerializer(
            job_type_list,
            many=True,
            context={
                'request': request
            }
        )

        return Response(serializer.data)
