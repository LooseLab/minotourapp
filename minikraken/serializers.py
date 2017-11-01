from rest_framework import serializers

from minikraken.models import MiniKraken
from minikraken.models import ParsedKraken


class MiniKrakenSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MiniKraken
        fields = ('run','read','krakenstatus','krakentaxid','krakenseqlen', 'krakenlca')


class ParsedKrakenSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = ParsedKraken
        fields = ('run','type_name','barcode_name','percentage','rootreads', 'directreads','rank','NCBItaxid', 'parent', 'sci_name', 'indentation', 'orderin')

