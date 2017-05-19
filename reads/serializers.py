from rest_framework import serializers

from reads.models import MinionRun, FastqRead, FastqReadType


class FastqReadTypeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqReadType
        fields = ('id', 'name',)
        read_only = ('id',)


class MinionRunSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinionRun
        fields = ('id', 'run_name', 'run_id', 'is_barcoded')
        read_only = ('id',)


class FastqReadSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqRead
        fields = ('id', 'read_id', 'read', 'channel', 'barcode', 'sequence', 'quality', 'is_pass',
                  'start_time')
        read_only = ('id',)
