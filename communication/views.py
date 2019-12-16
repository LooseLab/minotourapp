from django.shortcuts import render
from rest_framework.decorators import api_view
from rest_framework.response import Response

from communication.models import NotificationConditions
from communication.models import Message
from communication.serializers import MessageSerializer

# @api_view(['GET'])
# def new_messages_list(request):
#
#     queryset = Message.objects.filter(recipient=request.user).order_by('-created_date')
#     serializer = MessageSerializer(queryset, many=True, context={'request': request})
#     return Response(serializer.data)
#
#
# def message_details(request, pk):
#
#     message = Message.objects.get(pk=pk)
#
#     return render(
#         request,
#         'communication/message.html',
#         context={
#            'message': message,
#         }
#     )
from reads.models import Flowcell, Barcode
from reference.models import ReferenceInfo, ReferenceLine


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

    """
    print(kwargs)
    reference = kwargs.get("reference", None)
    chromosome = kwargs.get("chromosome", None)
    barcode = kwargs.get("barcode", None)
    set_repeat = kwargs.get("set_repeat", False)
    coverage_target = kwargs.get("cov_targ", 0)

    cond = NotificationConditions(notification_type=kwargs["condition"],
                                  flowcell=flowcell,
                                  creating_user=user,
                                  repeat=set_repeat,
                                  coverage_target=coverage_target,
                                  conditional_reference=reference,
                                  chromosome=chromosome,
                                  barcode=barcode,
                                  )
    return cond


@api_view(["GET", "POST"])
def get_or_create_conditions(request):
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
    condition_lookup = {"coverage": "cov", "newTarget": "targ", "mux": "mux"}

    if request.method == "POST":
        print(request.data)
        print(type(request.user))


        flowcell_pk = request.data.get("flowcell", None)

        flowcell = Flowcell.objects.get(pk=flowcell_pk)

        conditions = request.data.get("conditions", None)

        coverage_settings = request.data.get("coverage_sets", {})

        if flowcell is None:
            return Response("No flowcell submitted", status=400)
        if conditions is None:
            return Response("No conditions submitted", status=400)
        # Create conditions for each contained in the request body
        for condition, value in conditions.items():

            print(f"Condtion is {condition}")

            cov_targ = value if condition is "coverage" else None

            repeat = True

            condition_name = condition_lookup[condition]

            print(bool(coverage_settings))
            print(coverage_settings)
            print(condition)

            if condition == "coverage" and bool(coverage_settings):
                print("condition is coverage 2")

                repeat = False

                reference = ReferenceInfo.objects.get(pk=coverage_settings["reference"])
                chromosome = ReferenceLine.objects.get(pk=coverage_settings["chromosome"])
                kwargs = {"reference": reference, "chromosome": chromosome, "condition": condition_name,
                          "set_repeat": repeat, "cov_targ": cov_targ}
                print(f"barcodes are {coverage_settings['barcodes']}")

                if isinstance(coverage_settings["barcodes"], list):
                    print("making the conditions")
                    for barcode in coverage_settings["barcodes"]:
                        kwargs["barcode"] = Barcode.objects.get(pk=barcode)
                        cond = create_notification_conditions(flowcell, request.user, **kwargs)
                        print(cond)
                        cond.save()

            else:
                print("mux or target")
                kwargs = {"condition": condition_name, "set_repeat": True}
                cond = create_notification_conditions(flowcell, request.user, **kwargs)
                cond.save()
                print(cond)

        return Response("Successfully saved notification conditions.", status=200)


@api_view(["GET"])
def retrieve_messages(request):
    """
    Retrieve all the messages that have been sent to a user for display.
    Parameters
    ----------
    request: rest_framework.request.Request

    Returns
    -------

    """
    pass
