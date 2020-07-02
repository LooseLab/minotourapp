from collections import defaultdict

import pandas as pd
from django.db.models import F
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from alignment.models import PafSummaryCov
from communication.models import Message
from communication.models import NotificationConditions
from communication.serializers import MessageSerializer, NotificationSerialiser
from reads.models import JobMaster
from .forms import NotificationConditionForm


@api_view(["POST"])
def simple_condition_create(request):
    """
    Simple create using form
    Parameters
    ----------
    request

    Returns
    -------

    """
    print(request.data)
    form = NotificationConditionForm(request.data)
    if form.is_valid():
        new_condition = form.save(commit=False)
        new_condition.creating_user = request.user
        new_condition.save()

    return Response(status=status.HTTP_200_OK)


@api_view(["GET"])
def get_extant_condition_types(request, pk):
    """
    Get the extant condition types to disable the selection tabs
    Parameters
    ----------
    request: rest_framework.request.Request
        The request object in the body
    pk: int
        The primary key of the flowcell record in the database
    Returns
    -------
    list of str
        List of distinct string names of conditions that exist
    """
    results = (
        NotificationConditions.objects.filter(flowcell_id=pk)
        .values_list("notification_type", flat=True)
        .distinct()
    )
    return Response(results, status=status.HTTP_200_OK)


@api_view(["GET"])
def new_messages_list(request):
    """
    Get all the messages that have been sent by the twitter API
    Parameters
    ----------
    request: rest_framework.request.Request
        Django request framework request object.

    Returns
    -------

    """

    queryset = Message.objects.filter(recipient=request.user).order_by("-created_date")
    serializer = MessageSerializer(queryset, many=True, context={"request": request})
    return Response(serializer.data)


@api_view(["GET"])
def get_references_for_condtions(request, pk):
    """
    Get references for the coverage condition dropdowns
    Parameters
    ----------
    request: rest_framework.request.Request
        The request object for this AJAX request
    pk: int
        The flowcell record primary key
    Returns
    -------

    """
    b = defaultdict(list)
    for jobs in list(
        JobMaster.objects.filter(job_type_id=4).values(
            "reference_id",
            reference_name=F("reference__name"),
            contig_name=F("reference__referencelines__line_name"),
            contig_id=F("reference__referencelines__id"),
        )
    ):
        b[jobs["reference_name"]].append(
            (jobs["contig_name"], jobs["contig_id"], jobs["reference_id"])
        )

    return Response(b, status=status.HTTP_200_OK)


# def message_details(request, pk):
#     """
#     Returns the htm
#     Parameters
#     ----------
#     request
#     pk
#
#     Returns
#     -------
#
#     """
#     message = Message.objects.get(pk=pk)
#
#     return render(request, "communication/message.html", context={"message": message,})


def create_notification_conditions(flowcell, user, **kwargs):
    """
    Create and save a notification conditions
    Parameters
    ----------
    flowcell: reads.models.Flowcell
        The flowcell these notifications are being created for
    user:
        The user that we are creating this flowcell for
    kwargs: dict
        Key word arguments to create the flowcell with

    Returns
    -------
    cond: str
        Any error message (If any.)

    """
    print(kwargs)
    # reference = kwargs.get("reference", None)
    # chromosome = kwargs.get("chromosome", None)
    # barcode = kwargs.get("barcode", None)
    # set_repeat = kwargs.get("set_repeat", False)
    # coverage_target = kwargs.get("coverage_target", 0)

    lookup_dict = dict(NotificationConditions.NOTIFICATION_CHOICES)

    # Lookup the human readable notification type
    noti_type = lookup_dict[kwargs.get("notification_type")]

    # Check for my uniqueness

    check = NotificationConditions.objects.filter(
        flowcell=flowcell, creating_user=user, **kwargs
    )
    if check:
        return (
            f"An identical Condition already exists for the {noti_type}"
            f" condition you are trying to create."
        )

    # Now we check we don't have more than 2 non unique notifications of the same type, as long as they aren't Coverage
    # As many coverages as the user wants are allowed
    check_not_excess = 0
    if not kwargs.get("notification_type") == "cov":
        check_not_excess = NotificationConditions.objects.filter(
            flowcell=flowcell,
            creating_user=user,
            notification_type=kwargs.get("notification_type"),
        ).count()

    if check_not_excess > 2:

        return (
            f"Too many conditions of type {noti_type} already exist on this flowcell."
            f" Please delete an existing one."
        )

    if not check or not check_not_excess > 2:
        cond = NotificationConditions(flowcell=flowcell, creating_user=user, **kwargs)
        cond.save()


@api_view(["GET", "POST", "DELETE"])
def notification_conditions_list(request):
    """
    Get all the condition for a specific flowcell that a User has created
    Create the conditions that the user will be notified about if they are met. Chosen from notifications-manager.html.
    Parameters
    ----------
    request: rest_framework.request.Request
        The django rest framework request object. Contains the flowcell PK and a list of the conditions.

    Returns
    -------
    rest_framework.response.Response
    """

    if request.method == "DELETE":
        notification_pk = request.GET.get("pk", -1)
        try:
            deleted = NotificationConditions.objects.get(
                pk=int(notification_pk)
            ).delete()
            return Response(f"Deleted {deleted}", status=200)
        except Exception as e:
            return Response(e, status=500)
    elif request.method == "GET":
        flowcell_id = request.GET.get("flowcellId", -1)
        if flowcell_id == -1:
            return Response("Flowcell ID not provided.", status=400)
        queryset = NotificationConditions.objects.filter(
            creating_user=request.user, flowcell_id=flowcell_id
        )
        notifications_serialiser = NotificationSerialiser(queryset, many=True)
        return Response({"data": notifications_serialiser.data}, status=200)

@api_view(["GET"])
def get_coverage_summary(request, pk):
    """
    Get a list of all the chromosomes and their coverage.
    Parameters
    ----------
    request: rest_framework.request.Request
        Django rest framework object
    pk: int
        The primary key of the flowcell
    Returns
    -------
    If names is True - Dict of Reference to Contig, Contig to Barcodes, else return coverage per contig

    """
    # if the request is from the Dropdown for notifications
    if request.GET.get("names", False):
        queryset = JobMaster.objects.filter(
            flowcell__id=pk, job_type__name="minimap2"
        ).values_list(
            "reference__name",
            "reference__referencelines__line_name",
            "flowcell__runs__barcodes__name",
            "reference__id",
            "reference__referencelines__id",
            "flowcell__runs__barcodes__id",
        )

        if not queryset:
            return Response("No mapping tasks available on this flowcell.", status=404)

        reference_to_contig_dict = defaultdict(list)

        # create a dict of string keys and values to create the conditional drop downs
        # contig[0] - reference name
        # contig[1] - contig names
        # contig[2] - barcode names
        # contig[3] - reference pk
        # contig[4] - contig id
        # contig[5] - barcodes id
        for contig in queryset:

            contig_to_barcode_dict = defaultdict(list)

            reference_to_contig_dict[contig[0]] = [contig_to_barcode_dict, contig[3]]

            if not contig_to_barcode_dict[contig[1]]:

                contig_to_barcode_dict[contig[1]].append(contig[4])

            contig_to_barcode_dict[contig[1]].append((contig[2], contig[5]))

        return Response(reference_to_contig_dict, status=200)

    queryset = PafSummaryCov.objects.filter(job_master__flowcell__id=pk)

    if not queryset:
        return Response(
            "No alignment results found for this flowcell",
            status=status.HTTP_204_NO_CONTENT,
        )

    df = pd.DataFrame.from_records(
        queryset.values(
            "barcode_name",
            "chromosome__line_name",
            "total_yield",
            "reference_line_length",
        )
    )

    df["coverage"] = (
        df["total_yield"].div(df["reference_line_length"]).round(decimals=3)
    )

    df = df.drop(columns=["reference_line_length", "total_yield"])

    df.set_index("chromosome__line_name", inplace=True)

    dictdf = df.to_dict("index")

    result = {"data": dictdf}

    return Response(result, status=200)
