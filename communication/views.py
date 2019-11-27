from django.shortcuts import render
from rest_framework.decorators import api_view
from rest_framework.response import Response

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


@api_view(["GET", "POST"])
def get_or_create_conditions(request):
    """
    Get all the condition for a specific flowcell that a User has created
    Create the conditions that the user will be notified about if they are met. Chosen from notifications-manager.html.
    Parameters
    ----------
    request: django.core.handlers.wsgi.WSGIRequest

    Returns
    -------

    """
    if request.method == "POST":
        print(request.data)

        flowcell = request.data.get("flowcell", None)

        conditions = request.data.get("conditions", None)

        if flowcell is None:
            return Response("No flowcell submitted", status=400)
        if conditions is None:
            return Response("No conditions submitted", status=400)

        for condition in conditions:
            pass



        return Response("Hello", status=200)


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
