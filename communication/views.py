from collections import defaultdict

import pandas as pd
from django.db.models import F
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from alignment.models import PafSummaryCov
from communication.models import Message, NotificationConditionsBarcode
from communication.models import NotificationConditions
from communication.serializers import MessageSerializer, NotificationSerialiser
from reads.models import JobMaster, Barcode
from .forms import NotificationConditionForm


def create_condition_barcodes(barcodes, condition):
    """
    Create the Notification condition barcodes for all barcodes selected for coverage type conditions
    Parameters
    ----------
    barcodes: list
        List of barcode primary keys to creat for each condition
    condition: NotificationConditions
        The notification condition that we are creating

    Returns
    -------
    None
    """
    for barcode_pk in barcodes:
        barcode = Barcode.objects.get(pk=barcode_pk)
        NotificationConditionsBarcode.objects.get_or_create(
            barcode=barcode, barcode_name=barcode.name, condition=condition
        )


def save_condition(data, form, user):
    """
    Save the notification condition
    Parameters
    ----------
    data: dict
        The data that ahs been sent from post request
    form: communication.forms.NotificationConditionForm
        The form for this form
    user:
        django user
    Returns
    -------

    """
    new_condition = form.save(commit=False)
    new_condition.creating_user = user
    new_condition.save()
    if new_condition.notification_type == "cov":
        barcodes = list(map(int, data.get("barcodes", None).split(",")))
        create_condition_barcodes(barcodes, new_condition)


@api_view(["POST"])
def simple_condition_create(request):
    """
    Simple create using form
    Parameters
    ----------
    request: rest_framework.request.Request
        Post request from the communications page
    Returns
    -------

    """
    form = NotificationConditionForm(request.data)
    if form.is_valid():
        save_condition(request.data, form, request.user)
    else:
        # If we are doing genome level coverage, we have a chromosome of 0
        if form.data.get("chromosome", False) == "0":
            data = form.data.dict()
            data.pop("chromosome")
            form = NotificationConditionForm(data)
            if form.is_valid():
                save_condition(data, form, request.user)
            else:
                print(form.errors)
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
reads_
    Returns
    -------

    """
    query_columns = [
        "title",
        "created_date",
        "sender_first_name",
        "flowcell_name",
    ]
    draw = int(request.GET.get("draw", 0))
    search_value = request.GET.get("search[value]", "")
    start = int(request.GET.get("start", 0))
    length = int(request.GET.get("length", 10))
    end = start + length
    order_column = query_columns[int(request.GET.get("order[0][column]", ""))]
    # ascending descending
    order_dir = "-" if request.GET.get("order[0][dir]", "") == "desc" else ""
    queryset = Message.objects.filter(recipient=request.user)
    records_total = queryset.count()
    if search_value:
        queryset = queryset.filter(title__contains=search_value)
    queryset = queryset.order_by(f"{order_dir}{order_column}")[start:end]
    serializer = MessageSerializer(queryset, many=True, context={"request": request})
    result = {
        "draw": draw,
        "recordsTotal": records_total,
        "recordsFiltered": records_total,
        "data": serializer.data,
    }
    return Response(result)


@api_view(["GET"])
def get_references_for_conditions(request, pk):
    """
    Get references, contigs and barcodes for the coverage condition dropdowns
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
    barcodes = set()
    for contig in list(
        JobMaster.objects.filter(job_type_id=4, flowcell_id=pk)
        .values(
            "reference_id",
            "id",
            reference_name=F("reference__name"),
            contig_name=F("reference__reference_lines__line_name"),
            contig_id=F("reference__reference_lines__id"),
        )
        .distinct()
    ):
        b[contig["reference_name"]].append(
            (contig["contig_name"], contig["contig_id"], contig["reference_id"])
        )
        barcodes.update(
            PafSummaryCov.objects.filter(job_master_id=contig["id"])
            .exclude(barcode_name="No_barcode")
            .values_list("barcode_id", "barcode__name", "reference_pk")
        )
        b[f"{contig['reference_name']}_barcodes"] = list(barcodes)
    return Response(b, status=status.HTTP_200_OK)


@api_view(["GET", "DELETE"])
def notification_conditions_list(request):
    """
    Get all the condition for a specific flowcell that a User has created
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
            deleted = NotificationConditions.objects.get(pk=int(notification_pk))
            type = deleted.notification_type
            deleted.delete()
            return Response(type, status=200)
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
