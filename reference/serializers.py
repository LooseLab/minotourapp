from rest_framework import serializers

from reference.models import ReferenceInfo


class ReferenceInfoSerializer(serializers.ModelSerializer):
    """
    Serialiser for ReferenceInfos
    """
    class Meta:
        model = ReferenceInfo
        fields = (
            'id',
            'name',
            'length',
            'private',
        )
        read_only = (
            'id',
        )
