from reads.models import MinionRun

from rest_framework import serializers


class MinionRunSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinionRun
        fields = ('id', 'run_name', 'run_id', 'is_barcoded')
        read_only = ('id',)
