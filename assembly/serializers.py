from rest_framework import serializers

from assembly.models import GfaStore

class GfaStoreSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = GfaStore
        fields = (
            'id',
            'run',
            #'barcode',
            'nreads',
            'ncontigs',
            'maxlen',
            'totlen',
            'n50len'
        )
        read_only = (
            'id',
        )
