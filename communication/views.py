from django.shortcuts import render
from rest_framework.decorators import api_view
from rest_framework.response import Response

from communication.models import Message
from communication.serializers import MessageSerializer


@api_view(['GET'])
def new_messages_list(request):

    if request.method == 'GET':
        queryset = Message.objects.filter(recipient=request.user).order_by('-created_date')
        serializer = MessageSerializer(queryset, many=True, context={'request': request})
        return Response(serializer.data)


def message_details(request, pk):

    message = Message.objects.get(pk=pk)

    return render(
        request,
        'communication/message.html',
        context={
           'message': message,
        }
    )
