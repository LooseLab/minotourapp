from django.shortcuts import render
from rest_framework.decorators import api_view
from rest_framework.response import Response

from communication.models import NotificationConditions
from communication.models import Message
from communication.serializers import MessageSerializer, NotificationSerialiser
from reads.models import Flowcell, Barcode, UserOptions
from reference.models import ReferenceInfo, ReferenceLine


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


def message_details(request, pk):
    """
    Returns the htm
    Parameters
    ----------
    request
    pk

    Returns
    -------

    """
    message = Message.objects.get(pk=pk)

    return render(request, "communication/message.html", context={"message": message,})


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
        return f"An identical Condition already exists for the {noti_type}" \
               f" condition you are trying to create."

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
def get_create_delete_conditions(request):
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
    condition_lookup = {
        "coverage": "cov",
        "newTarget": "targ",
        "Mux alerts": "mux",
        "Voltage": "volt",
        "Temperature": "temp",
        "Warnings/Errors": "w/e",
        "Occupancy": "occu",
        "MinoTour": "mino",
    }

    if request.method == "POST":

        flowcell_pk = request.data.get("flowcell", None)

        flowcell = Flowcell.objects.get(pk=flowcell_pk)

        conditions = request.data.get("conditions", None)

        coverage_settings = request.data.get("coverage_sets", {})

        auto_choices = request.data.get("autoChoices", [])

        auto_values = request.data.get("rangeValues", [])

        twitter_settings = UserOptions.objects.get(owner=request.user)

        tweet = twitter_settings.tweet

        handle = twitter_settings.twitterhandle

        errors_list = []

        if flowcell is None:
            return Response("No flowcell submitted", status=400)
        if conditions is None:
            return Response("No conditions submitted", status=400)
        if not tweet and bool(handle):
            return Response(
                "Please provide Twitter permission in the profile section, and set a handle.",
                status=400,
            )
        # Create conditions for each contained in the request body
        for condition, value in conditions.items():

            cov_targ = value if condition == "coverage" else None

            condition_name = condition_lookup[condition]

            if condition == "coverage" and bool(coverage_settings) and value:

                repeat = False

                reference = ReferenceInfo.objects.get(pk=coverage_settings["reference"])
                chromosome = ReferenceLine.objects.get(
                    pk=coverage_settings["chromosome"]
                )
                kwargs = {
                    "conditional_reference": reference,
                    "chromosome": chromosome,
                    "notification_type": condition_name,
                    "repeat": repeat,
                    "coverage_target": cov_targ,
                }

                if isinstance(coverage_settings["barcodes"], list):

                    for barcode in coverage_settings["barcodes"]:
                        kwargs["barcode"] = Barcode.objects.get(pk=barcode)
                        cond = create_notification_conditions(
                            flowcell, request.user, **kwargs
                        )
                        errors_list.append(cond)

            # Mux or new target
        for condition in auto_choices:
            upper_limit = None
            lower_limit = None
            defaults = {"Occupancy": [40, 100], "Voltage": [-180, -160]}
            if condition in ["Occupancy", "Voltage"]:
                upper_limit = auto_values.get(condition, {}).get(
                    "upper", defaults[condition][1]
                )
                lower_limit = auto_values.get(condition, {}).get(
                    "lower", defaults[condition][0]
                )
            condition_name = condition_lookup[condition]
            kwargs = {
                "notification_type": condition_name,
                "repeat": True,
                "upper_limit": upper_limit,
                "lower_limit": lower_limit,
            }
            cond = create_notification_conditions(flowcell, request.user, **kwargs)
            errors_list.append(cond)
        return_list = ["Successfully created notification conditions.", "Issues:"]
        # Remove None, as successfully create notification conditions function returns none.
        errors_list = list(filter(None, errors_list))
        # If we have no errors
        if not errors_list:
            return_list[1] = "Issues: None."
        return_list.extend(errors_list)
        return_list = list(filter(None, return_list))
        return_string = "<br/><br/>".join(return_list)
        return Response(return_string, status=200)

    elif request.method == "GET":

        flowcell_id = request.GET.get("flowcellId", -1)

        if flowcell_id == -1:
            return Response("Flowcell ID not provided.", status=400)

        queryset = NotificationConditions.objects.filter(
            creating_user=request.user, flowcell_id=flowcell_id
        )

        notifications_serialiser = NotificationSerialiser(queryset, many=True)

        return Response({"data": notifications_serialiser.data}, status=200)

    elif request.method == "DELETE":

        notification_pk = request.GET["ntpk"]

        try:

            deleted = NotificationConditions.objects.get(pk=int(notification_pk)).delete()
            return Response(f"Deleted {deleted}", status=200)

        except Exception as e:
            return Response(e, status=500)

