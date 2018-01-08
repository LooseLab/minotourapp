from rest_framework import serializers

from alignment.models import PafStore
from alignment.models import PafRoughCov
from alignment.models import PafSummaryCov
from reads.models import BarcodeGroup

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
            'barcodegroup',
            'reference',
            'ref_name',
            'ref_len',
            'chromosome',
            'chrom_name',
            'chrom_len',
            'read_count',
            'cumu_length',
            'chrom_cover',
            'avg_read_len',
            'read_type_name',
            'barcode_name',
            'barcode_group_name'
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