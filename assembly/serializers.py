from rest_framework import serializers

from assembly.models import GfaSummary

class GfaSummarySerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = GfaSummary
        fields = (
            'id',
            'run',
            'barcode_name',
            'type_name',
            'nreads',
            'ncontigs',
            'maxlen',
            'totlen',
            'n50len',
            'timecreated'
        )
        read_only = (
            'id',
        )
