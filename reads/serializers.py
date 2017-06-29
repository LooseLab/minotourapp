from rest_framework import serializers

from reads.models import MinIONRun, \
    FastqRead, \
    FastqReadType, \
    MinION, \
    MinIONEvent, \
    MinIONEventType, \
    MinIONScripts, \
    MinIONStatus,   \
    MinIONRunStatus, \
    MinIONRunStats, \
    MinIONmessages, \
    MinIONControl


class MinIONRunSerializer(serializers.HyperlinkedModelSerializer):
    # reads = serializers.HyperlinkedRelatedField(many=True,
    #                                            view_name='fastqread-detail',
    #                                            read_only=True)

    class Meta:
        model = MinIONRun
        # fields = ('url', 'run_name', 'run_id', 'is_barcoded', 'reads')
        fields = ('url', 'sample_name', 'minKNOW_version', 'minKNOW_flow_cell_id', 'run_name', 'run_id', 'is_barcoded','minION')#, 'barcodes')
        read_only = ('id','sample_name','minKNOW_version', 'minKNOW_flow_cell_id')#, 'barcodes')

class MinIONRunStatsSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinIONRunStats
        fields = (
            #'__all__'
            'id',
            'minION',
            'run_id',
            'sample_time',
            'event_yield',
            'asic_temp',
            'heat_sink_temp',
            'voltage_value',
            'mean_ratio',
            'open_pore',
            'in_strand',
            'multiple',
            'unavailable',
            'unknown',
            'adapter',
            'pending_mux_change',
            'unclassified',
            'below',
            'unblocking',
            'above',
            'good_single',
            'saturated',
            'inrange',
            'strand',
            'occupancy',
            'minKNOW_read_count',
            'minKNOW_histogram_values',
            'minKNOW_histogram_bin_width',
        )
        read_only = ('id','occupancy',)

class MinIONSerializer(serializers.HyperlinkedModelSerializer):
    #currentrunstats = MinIONRunStatsSerializer(many=True, read_only=True)

    class Meta:
        model = MinION
        fields = ('url', 'space_available','minKNOW_version', 'minION_name', 'status', 'computer','sample_name','run_status','flow_cell_id','run_name','total_drive_space','space_till_shutdown','warnings','last_run','currentscript','event_yield','voltage_value',)
        read_only = ('status','space_available','minKNOW_version','computer','sample_name','run_status','flow_cell_id','run_name','total_drive_space','space_till_shutdown','warnings','last_run','currentscript','event_yield','voltage_value',)

class FastqReadTypeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqReadType
        fields = ('url', 'name',)
        read_only = ('id',)

class MinIONmessagesSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinIONmessages
        fields = ('minION', 'run_id', 'minKNOW_message', 'minKNOW_severity', 'minKNOW_message_timestamp', 'minKNOW_identifier',)
        read_only = ('id',)


class MinIONControlSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinIONControl
        fields = ('id', 'minION', 'job', 'custom','complete')
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


class MinIONScriptsSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinIONScripts
        fields = ('id', 'minION', 'identifier', 'name', 'experiment_type', 'base_calling', 'flow_cell', 'kit')
        read_only = ('id',)

class MinIONStatusSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinIONStatus
        fields = (
            'id',
            'minION',
            'minKNOW_status',
            'minKNOW_current_script',
            'minKNOW_sample_name',
            'minKNOW_exp_script_purpose',
            'minKNOW_flow_cell_id',
            'minKNOW_run_name',
            'minKNOW_hash_run_id',
            'minKNOW_script_run_id',
            'minKNOW_real_sample_rate',
            'minKNOW_asic_id',
            'minKNOW_total_drive_space',
            'minKNOW_disk_space_till_shutdown',
            'minKNOW_disk_available',
            'minKNOW_warnings',

        )
        read_only = ('id',)




class MinIONRunStatusSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinIONRunStatus
        fields = (
            'id',
            'minION',
            'minION_name',
            #'minKNOW_status',
            'minKNOW_current_script',
            'minKNOW_sample_name',
            'minKNOW_exp_script_purpose',
            'minKNOW_flow_cell_id',
            'minKNOW_run_name',
            'minKNOW_version',
            'run_id',
            'minKNOW_hash_run_id',
            'minKNOW_script_run_id',
            'minKNOW_real_sample_rate',
            'minKNOW_asic_id',
            'minKNOW_start_time',
            'minKNOW_colours_string'

            #'minKNOW_total_drive_space',
            #'minKNOW_disk_space_till_shutdown',
            #'minKNOW_warnings',
        )
        read_only = ('id','minION_name')


