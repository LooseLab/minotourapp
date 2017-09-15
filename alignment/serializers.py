from rest_framework import serializers

from alignment.models import PafStore
from alignment.models import PafRoughCov


class PafStoreSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = PafStore
        fields = '__all__'


class PafRoughCovSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = PafRoughCov
        fields = (
            'id',
            'run',
            'barcode',
            'reference',
            'chromosome',
            'position',
            'incdel',
        )
        read_only = (
            'id',
        )