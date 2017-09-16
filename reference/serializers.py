from rest_framework import serializers

from reference.models import ReferenceInfo
from reference.models import ReferenceLine

class ReferenceInfoSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = ReferenceInfo
        fields = (
            #'id',
            'reference_name',
            'filename',
            'bwa_index_file_location',
            'minimap2_index_file_location',
            'totalrefleN',
            'private',
            'users'
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