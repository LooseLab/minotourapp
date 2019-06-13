from django.http import JsonResponse, HttpResponse
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.parsers import JSONParser
from rest_framework.response import Response

from jobs.models import JobType, JobMaster
from jobs.serializers import JobMasterSerializer, JobMasterInsertSerializer
from reads.models import Run
from reference.models import ReferenceInfo
from web.forms import TaskForm
from reads.models import Flowcell


@api_view(["GET", "POST"])
def task_list(request):
    """
    return aLL the tasks active under one flowcell
    :param request:
    :return:
    """
    print(request.GET)
    if request.method == 'GET':

        search_criteria = request.GET.get('search_criteria', 'flowcell')

        if search_criteria == 'flowcell':
            flowcell_id = request.GET.get("search_value", "")  # TODO have to test for empty search value
            task_list = JobMaster.objects.filter(flowcell__id=int(flowcell_id)).exclude(job_type__name="Other")
            serializer = JobMasterSerializer(task_list, many=True)
            result = {
                "data": serializer.data
            }
            return JsonResponse(result)

        else:
            return JsonResponse({"data": []})

    else:
        # If it's a metagenomics tasks, and it's from the web app not the client, so target set isn't in the dict keys
        if request.data["job_type"] == "10" and "target_set" not in request.data.keys():
            # Set the target set id to the provided reference ID
            request.data["target_set"] = request.data["reference"]
            # Remove the reference ID from the request body
            request.data["reference"] = None

        if request.data["job_type"] == "11":
            request.data["reference"] = None

        # Check to see if we have a string as the flowcell, not an Int for a PK based lookup
        if type(request.data["flowcell"]) is str:
            try:
                flowcell = Flowcell.objects.get(name=request.data["flowcell"])
                request.data["flowcell"] = flowcell.id
            except Flowcell.DoesNotExist as e:
                try:
                    flowcell = Flowcell.objects.get(pk=int(request.data["flowcell"]))
                except Flowcell.DoesNotExist as e1:
                    print("Exception: {}".format(e))
                    return JsonResponse({'error_messages': "Exception: {}".format(e)}, status=500)

        if "reference" in request.data.keys() and type(request.data["reference"]) is str:
            try:
                reference = ReferenceInfo.objects.get(name=request.data["reference"])
                request.data["reference"] = reference.id
            except ReferenceInfo.DoesNotExist as e:
                try:
                    reference = ReferenceInfo.objects.get(pk=int(request.data["reference"]))
                except ReferenceInfo.DoesNotExist as e1:
                    print("Exception: {}".format(e))
                    return JsonResponse({'error_messages': "Exception: {}".format(e)}, status=500)

        serializer = JobMasterInsertSerializer(data=request.data)
        if serializer.is_valid():
            task = serializer.save()

            response_data = {}
            response_data['message'] = "Create task successful!"
            response_data['pk'] = task.id

            return JsonResponse(response_data)

        else:

            return JsonResponse({'error_messages': serializer.errors}, status=500)


@api_view(['GET'])
def task_types_list(request):

    if request.GET.get("cli", False):
        tasks = ["Metagenomics", "Assembly", "Minimap2"]
        queryset = JobType.objects.filter(name__in=tasks)
    else:
        queryset = JobType.objects.filter(private=False)

    result = []

    for record in queryset:
        task = {
            'id': record.id,
            'name': record.name,
            'description': record.description
        }

        result.append(task)

    return JsonResponse({
        'data': result
    })


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

