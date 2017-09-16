from rest_framework import serializers

from alignment.models import PafStore
from alignment.models import PafRoughCov

from reference.models import ReferenceLine
from reference.models import ReferenceInfo

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
            'p',
            'i',
        )
        read_only = (
            'id',
        )

class PafRoughCovChromSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = PafRoughCov
        fields = (
            #'id',
            #'run',
            #'barcode',
            #'reference',
            #'chromosome',
            'p',
            'i',
        )
        read_only = (
            'id',
        )