from rest_framework import serializers

from alignment.models import PafRoughCov, PafStore


class PafStoreSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = PafStore
        fields = '__all__'


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


class PafRoughCovSerializer(serializers.ModelSerializer):

    class Meta:

        model = PafRoughCov

        fields = (
            'p',
            'i',
        )
