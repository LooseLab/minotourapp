from django.http import JsonResponse
from rest_framework.decorators import api_view
from rest_framework.response import Response

from jobs.models import JobType, JobMaster
from jobs.serializers import JobMasterSerializer
from reads.serializers import JobTypeSerializer


@api_view(["GET"])
def task_list(request):
    """
    return aLL the tasks active under one flowcell
    :param request:
    :return:
    """
    search_criteria = request.GET.get('search_criteria', 'flowcell')

    if search_criteria == 'flowcell':
        flowcell_id = request.GET.get("search_value", "") # TODO have to test for empty search value
        task_list = JobMaster.objects.filter(flowcell__id=int(flowcell_id))
        serializer = JobMasterSerializer(task_list, many=True)
        result = {
            "data": serializer.data
        }
        return JsonResponse(result)
    else:
        return JsonResponse({"data": []})

@api_view(['GET'])
def task_types_list(request):
    if request.method == 'GET':
        queryset = JobType.objects.filter(private=False)
        serializer = JobTypeSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)

@api_view(['POST'])
def set_task_detail_all(request, pk):
    """We need to check if a job type already exists - if so we are not going to add another."""
    if request.method == 'POST':
        jobtype = JobType.objects.get(name=request.data["job"])
        print(jobtype)
        reference = ""
        if request.data["reference"] != "null":
            reference = ReferenceInfo.objects.get(reference_name=request.data["reference"])
            print(reference)
        minionrun = Run.objects.get(id=pk)
        print(minionrun)
        print(request.data)
        print(jobtype, reference, minionrun)
        jobmasters = JobMaster.objects.filter(run=minionrun).filter(job_type=jobtype)
        print("Jobmasters", jobmasters)
        if len(jobmasters) > 0:
            # return Response("Duplicate Job attempted. Not allowed.", status=status.HTTP_400_BAD_REQUEST)
            return Response("Duplicate Job attempted. Not allowed.", status=status.HTTP_200_OK)
        else:
            newjob = JobMaster(run=minionrun, job_type=jobtype, last_read=0, read_count=0, complete=False,
                               running=False)

            print("trying to make a job", newjob)

            if request.data["reference"] != "null":
                # if len(reference)>0:
                newjob.reference = reference
            try:
                newjob.save()
                print("job created")
            except Exception as e:
                print(e)
                print("error")
            return Response("Job Created.", status=status.HTTP_200_OK)

@api_view(['GET'])
def tasks_detail_all(request, pk):

    queryset = JobType.objects.filter(private=False).filter(type_name__type_name__in=['run', ])
    minionrun = Run.objects.get(pk=pk)

    result = []
    # print (queryset)
    for jobtype in queryset:
        # print (jobtype.type_name.all())
        obj = {}
        obj.update({
            'name': jobtype.name,
            'description': jobtype.description,
            'long_description': jobtype.long_description,
            'read_count': jobtype.readcount,
            'reference': jobtype.reference,
            'transcriptome': jobtype.transcriptome
        })

        jobmasterlist = JobMaster.objects.filter(run=minionrun).filter(job_type=jobtype)

        if len(jobmasterlist) > 0:
            obj2 = {}
            if jobmasterlist[0].reference:
                reference_name = jobmasterlist[0].reference.reference_name
            else:
                reference_name = ''
            obj2.update({
                'reference': reference_name,
                'last_read': jobmasterlist[0].last_read,
                'read_count': jobmasterlist[0].read_count,
                'temp_file': jobmasterlist[0].tempfile_name,
                'complete': jobmasterlist[0].complete,
                'running': jobmasterlist[0].running
            })

            obj.update({
                'job_details': obj2
            })

        result.append(obj)

    return HttpResponse(json.dumps(result), content_type="application/json")

