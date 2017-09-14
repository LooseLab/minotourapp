from rest_framework import serializers

from alignment.models import PafStore


class PafStoreSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = PafStore
        fields = ('__all__',)