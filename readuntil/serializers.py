from rest_framework import serializers

from readuntil.models import RejectedFastqRead


class RejectedFastqReadSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:

        model = RejectedFastqRead

        fields = (
            'run',
            'read_id',
        )
