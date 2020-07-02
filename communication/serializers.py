from rest_framework import serializers

from communication.models import Message, NotificationConditions


class MessageSerializer(serializers.ModelSerializer):

    recipient_email = serializers.ReadOnlyField(source="recipient.email", default="Unknown")
    sender_email = serializers.ReadOnlyField(source="sender.email", default="Unknown")
    sender_first_name = serializers.ReadOnlyField(source="sender.first_name", default="Unknown")
    flowcell_name = serializers.ReadOnlyField(source="flowcell.name", default="Unknown")

    class Meta:

        model = Message
        fields = (
            'uuidstr',
            'recipient_email',
            'sender_email',
            'sender_first_name',
            'title',
            'is_read',
            'created_date',
            'flowcell',
            'run',
            'flowcell_name'
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
        #


class NotificationSerialiser(serializers.ModelSerializer):
    """
    Serialiser used to retrieve the Notifications to be displayed in a table.
    """
    ref_name = serializers.ReadOnlyField(source="conditional_reference.name", default="N/A")
    chrom_line_name = serializers.ReadOnlyField(source="chromosome.line_name", default="N/A")
    # Get human readable name of Choice
    notification_type = serializers.CharField(read_only=True, source='get_notification_type_display')
    class Meta:
        model = NotificationConditions
        fields = (
            "id",
            "notification_type",
            "completed",
            "coverage_target",
            "upper_limit",
            "lower_limit",
            "ref_name",
            "chrom_line_name"

        )
    def to_representation(self, instance):
        my_fields = {'upper_limit', 'lower_limit', 'coverage_target'}
        data = super().to_representation(instance)
        for field in my_fields:
            try:
                if not data[field]:
                    data[field] = "N/A"
            except KeyError:
                pass
        return data
