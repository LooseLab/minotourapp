from rest_framework import serializers

from minikraken.models import MiniKraken


class MiniKrakenSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MiniKraken
        fields = '__all__'