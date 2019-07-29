from django.http import JsonResponse, HttpResponse
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from jobs.models import JobType, JobMaster
from jobs.serializers import JobMasterSerializer, JobMasterInsertSerializer
from reads.models import Run
from reference.models import ReferenceInfo
from reads.models import Flowcell
from web.delete_tasks import delete_metagenomics_task, delete_alignment_task, delete_expected_benefit_task
import time

@api_view(["GET", "POST"])
def get_or_create_tasks(request):
    """
    API endpoint for dealing with fetching or creating new JobMasters
    :param request: The django rest framework request body - if a GET request should contain
    search_criteria and search_value params
    :return: A dictionary containing all the matching JobMasters on this flowcell
    """

    if request.method == 'GET':

        search_criteria = request.GET.get('search_criteria', 'flowcell')

        if search_criteria == 'flowcell':

            flowcell_id = request.GET.get("search_value", -1)
            # If there was no flowcell ID provided
            if flowcell_id == -1:
                return Response("No flowcell ID provided", status=status.HTTP_400_BAD_REQUEST)
            # Get all the JobMasters for this flowcell, excluding job type other - which are for background sub tasks
            tasks_list = JobMaster.objects.filter(flowcell__id=int(flowcell_id)).exclude(job_type__name="Other")
            # Serialise the data to a python object
            serializer = JobMasterSerializer(tasks_list, many=True)
            # Check the icon required on the pause control part of the table
            for data in serializer.data:
                if data["paused"]:
                    data["icon"] = "play"
                    data["iconText"] = "Play"
                else:
                    data["icon"] = "pause"
                    data["iconText"] = "Pause"

            result = {
                "data": serializer.data
            }

            return JsonResponse(result)

        else:

            return Response("Bad search criteria provided", status=status.HTTP_400_BAD_REQUEST)

    # Its a post request to create a new task
    else:
        # Check to see if we have a string as the flowcell, not an Int for a PK based lookup
        # This is for starting a job from the client

        # If the job isn't EB or minimap2
        if int(request.data["job_type"]) not in [4, 15]:

            request.data["reference"] = None

        if type(request.data["flowcell"]) is str:

            try:
                # Get the flowcell by name
                flowcell = Flowcell.objects.get(name=request.data["flowcell"])

                request.data["flowcell"] = flowcell.id
            # If that doesn't work
            except Flowcell.DoesNotExist as e:

                try:
                    # Get the flowcell by pk but change the string number to an int
                    flowcell = Flowcell.objects.get(pk=int(request.data["flowcell"]))

                    request.data["flowcell"] = flowcell.id

                except Flowcell.DoesNotExist:

                    print("Exception: {}".format(e))

                    return Response("Flowcell not found. Please contact server admin.", status=500)

        # For client side job starting, check if we have a reference and it's a string
        if "reference" in request.data.keys() and type(request.data["reference"]) is str:

            try:
                # If a reference hasn't been selected.
                if request.data["reference"] is "":
                    return Response("Reference does not exist - Please select a reference",
                                        status=500)

                # Try to get the reference by name
                reference = ReferenceInfo.objects.get(name=request.data["reference"])

                request.data["reference"] = reference.id

            except ReferenceInfo.DoesNotExist as e:

                try:

                    reference = ReferenceInfo.objects.get(pk=int(request.data["reference"]))

                    request.data["reference"] = reference.id

                except ReferenceInfo.DoesNotExist:

                    print("Exception: {}".format(e))

                    return Response("Reference not found Please contact server admin", status=500)
        # Serialise the data to a Django savable object
        print(request.data)
        serializer = JobMasterInsertSerializer(data=request.data)

        # If the serialiser is valid
        if serializer.is_valid():

            task = serializer.save()

            response_data = {}

            response_data['message'] = "Task created successfully!"

            response_data['pk'] = task.id

            return JsonResponse(response_data)

        else:

            return JsonResponse({'error_messages': serializer.errors}, status=500)


@api_view(['GET'])
def task_types_list(request):
    """
    Get the list of task types to display for Job creation options
    :param request: The request body object
    :return: A list of dictionaries containing the task id, name and description.
    """

    # If it's a request from the client
    if request.GET.get("cli", False):
        # These are the available tasks
        tasks = ["Metagenomics", "Assembly", "Minimap2"]
        # Get the tasks
        queryset = JobType.objects.filter(name__in=tasks)

    else:
        # Otherwise it's for the site, so make all that are public available
        queryset = JobType.objects.filter(private=False)

    # Create the result list
    result = [{
        'id': record.id,
        'name': record.name,
        'description': record.description
    } for record in queryset
    ]

    return JsonResponse({
        'data': result
    })


# @api_view(['POST'])
# def set_task_detail_all(request, pk):
#     # TODO ARE WE USING THIS AT ALL
#     """We need to check if a job type already exists - if so we are not going to add another."""
#     if request.method == 'POST':
#         jobtype = JobType.objects.get(name=request.data["job"])
#         print(jobtype)
#         reference = ""
#         if request.data["reference"] != "null":
#             reference = ReferenceInfo.objects.get(reference_name=request.data["reference"])
#             print(reference)
#         minionrun = Run.objects.get(id=pk)
#         print(minionrun)
#         print(request.data)
#         print(jobtype, reference, minionrun)
#         jobmasters = JobMaster.objects.filter(run=minionrun).filter(job_type=jobtype)
#         print("Jobmasters", jobmasters)
#         if len(jobmasters) > 0:
#             # return Response("Duplicate Job attempted. Not allowed.", status=status.HTTP_400_BAD_REQUEST)
#             return Response("Duplicate Job attempted. Not allowed.", status=status.HTTP_200_OK)
#         else:
#             newjob = JobMaster(run=minionrun, job_type=jobtype, last_read=0, read_count=0, complete=False,
#                                running=False)
#
#             print("trying to make a job", newjob)
#
#             if request.data["reference"] != "null":
#                 # if len(reference)>0:
#                 newjob.reference = reference
#             try:
#                 newjob.save()
#                 print("job created")
#             except Exception as e:
#                 print(e)
#                 print("error")
#             return Response("Job Created.", status=status.HTTP_200_OK)


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


@api_view(["POST"])
def task_control(request):
    """
    Endpoint for clearing a task and resetting it to start form the beginning
    :param request: Request object, contains the Id of the task keyed to flowcellJobId
    :return: A status reflecting the state of the executed code
    """
    # Get the task object from django ORM
    job_master = JobMaster.objects.get(pk=request.data["flowcellJobId"])

    action_type = request.data["actionType"]

    lookup_action_type = {1: "Reset", 2: "Pause", 3: "Delete"}

    action = lookup_action_type[action_type]

    unrecognised_action_message = "Apologies, but this action type was not recognised." \
                                  " It may have not been implemented yet."

    if job_master.job_type.name == "ChanCalc":
        if action == "Reset":

            running = job_master.running

            seconds_counter = 0

            while running:

                running = job_master.running
                time.sleep(1)
                seconds_counter += 1

                if seconds_counter is 25:
                    job_master.running = False
                    job_master.save()

            # reset the values, skips existing values
            job_master.read_count = 0
            job_master.last_read = 0
            job_master.save()
            # set a return message
            return_message = f"Successfully reset ChanCalc task, id: {job_master.id}"

        elif action == "Pause":

            if job_master.paused:
                job_master.paused = False
            if not job_master.paused:
                job_master.paused = True

            job_master.save()
            return_message = f"Successfully paused ChanCalc task, id: {job_master.id}"

        else:
            return Response(unrecognised_action_message, status=500)

    elif job_master.job_type.name == "UpdateFlowcellDetails":
        if action == "Reset":
            # reset the values, skips existing values
            running = job_master.running

            seconds_counter = 0

            while running:
                running = job_master.running
                time.sleep(1)
                seconds_counter += 1
                if seconds_counter is 30:
                    job_master.running = False
                    job_master.save()

            job_master.read_count = 0
            job_master.last_read = 0
            job_master.save()
            # set a return message
            return_message = f"Successfully reset UpdateFlowcellDetails task, id: {job_master.id}"

        elif action == "Pause":
            if job_master.paused:
                job_master.paused = False
            else:
                job_master.paused = True

            job_master.save()
            return_message = f"Successfully paused UpdateFlowcellDetails task, id: {job_master.id}"

        else:

            return Response(unrecognised_action_message, status=500)

    elif job_master.job_type.name == "Metagenomics":
        if action == "Reset":
            delete_metagenomics_task.delay(job_master.id, True)

            return_message = f"Successfully started to reset this metagenomics task, id: {job_master.id}." \
                             f" Clearing previous data may take a while, please be patient!"

        elif action == "Pause":
            if job_master.paused:
                job_master.paused = False
            else:
                job_master.paused = True

            job_master.save()
            return_message = f"Successfully paused metagenomics task reset, id: {job_master.id}."

        elif action == "Delete":
            delete_metagenomics_task.delay(job_master.id)
            return_message = f"Successfully started deletion of this metagenomics task, id: {job_master.id}." \
                             f" Clearing previous data may take a while, please be patient!"

        else:
            return Response(unrecognised_action_message, status=500)

    elif job_master.job_type.name == "Minimap2":
        if action == "Reset":
            delete_alignment_task.delay(job_master.id, True)
            return_message = f"Successfully began reset of minimap2 task, id: {job_master.id}" \
                             f" Clearing previous data may take a while, please be patient!"

        elif action == "Pause":
            if job_master.paused:
                job_master.paused = False
            else:
                job_master.paused = True
            job_master.save()
            return_message = f"Successfully paused minimap2 task, id: {job_master.id}"

        elif action == "Delete":
            delete_alignment_task.delay(job_master.id)
            return_message = f"Successfully began deletion of minimap2 task, id: {job_master.id}" \
                             f" Clearing previous data may take a while, please be patient!"
        else:
            return Response(unrecognised_action_message, status=500)

    elif job_master.job_type.name == "ExpectedBenefit":
        if action == "Reset":
            return_message = f"Successfully reset this Expected benefit task, id: {job_master.id}"

        elif action == "Pause":
            if job_master.paused:
                job_master.paused = False
            else:
                job_master.paused = True
            return_message = f"Successfully paused this Expected benefit task, id: {job_master.id}"

        elif action == "Delete":
            return_message = f"Successfully deleted this Expected benefit task, id: {job_master.id}"

        else:
            return Response(unrecognised_action_message, status=500)

    else:
        raise NotImplementedError

    return Response(return_message, status=200)

