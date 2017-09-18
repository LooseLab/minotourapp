from rest_framework import serializers

from alignment.models import PafStore
from alignment.models import PafRoughCov
from alignment.models import PafSummaryCov

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
            #'p',
            #'i',
        )
        read_only = (
            'id',
        )

class PafSummaryCovSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = PafSummaryCov
        fields = (
            'id',
            'run',
            'read_type',
            'barcode',
            'reference',
            'chromosome',
            'read_count',
            'cumu_length',
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

class PafRoughCovChromSerializerCount(serializers.HyperlinkedModelSerializer):
    sumcount = serializers.IntegerField()
    reference = serializers.IntegerField()

    class Meta:
        model = PafRoughCov
        fields = (
            #'id',
            #'run',
            #'barcode',
            'reference',
            #'chromosome',
            #'p',
            #'i',
            'sumcount',

        )
        read_only = (
            'id',
        )