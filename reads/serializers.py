from rest_framework import serializers

from reads.models import MinIONRun, FastqRead, FastqReadType, MinION, MinIONEvent, MinIONEventType


class MinIONRunSerializer(serializers.HyperlinkedModelSerializer):
    # reads = serializers.HyperlinkedRelatedField(many=True,
    #                                            view_name='fastqread-detail',
    #                                            read_only=True)

    class Meta:
        model = MinIONRun
        # fields = ('url', 'run_name', 'run_id', 'is_barcoded', 'reads')
        fields = ('url', 'run_name', 'run_id', 'is_barcoded')


class MinIONSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinION
        fields = ('url', 'minION_name',)

class FastqReadTypeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqReadType
        fields = ('url', 'name',)
        read_only = ('id',)

class MinIONEventTypeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinIONEventType
        fields = ('url', 'name',)
        read_only = ('id',)

class FastqReadSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqRead
        fields = ('url', 'read_id', 'read', 'channel', 'barcode', 'sequence', 'quality', 'is_pass',
                  'start_time', 'run_id', 'type')


class FastqReadNameSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqRead
        fields = ('read_id',)





class MinIONEventSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinIONEvent
        fields = ('id', 'computer_name','minION', 'event', 'datetime')
        read_only = ('id',)