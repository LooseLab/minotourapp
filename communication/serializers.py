from rest_framework import serializers

from communication.models import Message


class MessageSerializer(serializers.ModelSerializer):

    recipient_email = serializers.ReadOnlyField(source="recipient.email")
    sender_email = serializers.ReadOnlyField(source="sender.email")

    class Meta:

        model = Message
        fields = (
            'uuidstr',
            'recipient_email',
            'sender_email',
            'title',
            'is_read',
            'created_date',
        )


class MessageSerializerComplete(serializers.ModelSerializer):

    recipient_email = serializers.ReadOnlyField(source="recipient.email")
    sender_email = serializers.ReadOnlyField(source="sender.email")

    class Meta:

        model = Message
        fields = (
            'id',
            'recipient_email',
            'sender_email',
            'title',
            'content',
        )
