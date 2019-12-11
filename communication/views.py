from django.shortcuts import render
from rest_framework.decorators import api_view
from rest_framework.response import Response

from communication.models import NotificationConditions, ConditionalChromosomes, ConditionalBarcodes
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
        print(type(request))

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

            cov_targ = value if condition is "coverage" else None

            repeat = True

            reference = None

            if condition is "coverage" and coverage_settings:
                repeat = False
                reference = ReferenceInfo.objects.get(name=coverage_settings["reference"])

            cond = NotificationConditions(notification_type=condition_lookup[condition],
                                          flowcell=flowcell,
                                          creating_user=request.user,
                                          repeat=repeat,
                                          coverage_target=cov_targ,
                                          reference_file=reference
                                          )
            cond.save()
        if "coverage" in conditions:
            pass
            # a = {"barcodes": (ConditionalBarcodes, Barcode, ), "chromosome": (ConditionalChromosomes, ReferenceLine)}
            # for name, settings in coverage_settings.items():
            #     if isinstance(settings, list):
            #         for itemy in settings:
            #             foreign = a[name][1].objects.get(re)
            #             a[name][0](
            #
            #             )



        return Response("Successfully saved notification conditions.", status=200)


@api_view(["GET"])
def retrieve_messages(request):
    """
    Retrieve all the messages that have been sent to a user for display.
    Parameters
    ----------
    request: django.core.handlers.wsgi.WSGIRequest

    Returns
    -------

    """
    pass
