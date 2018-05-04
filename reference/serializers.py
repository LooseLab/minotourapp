from rest_framework import serializers

from reference.models import ReferenceInfo, ReferenceLine


class ReferenceInfoSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:

        model = ReferenceInfo

        fields = (
            'id',
            'url',
            'reference_name',
            'totalrefleN',
            'private',
        )
        read_only = (
            'id',
        )


class ReferenceLineSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = ReferenceLine
        fields = (
            'id',
            'reference',
            'line_name',
            'chromosome_length',
        )
        read_only = (
            'id',
        )
