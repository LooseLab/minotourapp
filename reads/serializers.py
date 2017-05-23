from rest_framework import serializers

from reads.models import MinionRun, FastqRead, FastqReadType


class MinionRunSerializer(serializers.HyperlinkedModelSerializer):
    # reads = serializers.HyperlinkedRelatedField(many=True,
    #                                            view_name='fastqread-detail',
    #                                            read_only=True)

    class Meta:
        model = MinionRun
        # fields = ('url', 'run_name', 'run_id', 'is_barcoded', 'reads')
        fields = ('url', 'run_name', 'run_id', 'is_barcoded')


class FastqReadTypeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqReadType
        fields = ('url', 'name',)
        read_only = ('id',)


class FastqReadSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqRead
        fields = ('url', 'read_id', 'read', 'channel', 'barcode', 'sequence', 'quality', 'is_pass',
                  'start_time', 'run_id', 'type')
