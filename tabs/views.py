from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response
from reads.models import MinIONRun
#from tabs.models import RunID
#from .serializers import RunIDSerializer

'''
@api_view(['GET'])
def tabs(request):
    """
    List of all runs with associated tabs
    """
    if request.method == 'GET':
        queryset = RunID.objects.all()
        serializer = RunIDSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)


@api_view(['GET'])
def tabs_details(request, pk):
    """
    Retrieve tabs associated with a runID.
    """
    return Response(RunID.objects.filter(run_id=pk).values_list('tab_id__name', flat=True))
'''