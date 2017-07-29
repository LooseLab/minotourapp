from rest_framework import serializers

from reads.models import FastqRead, RunSummaryBarCode, RunStatistic, RunStatisticBarcode
from reads.models import FastqReadType
from reads.models import MinION
from reads.models import MinIONControl
from reads.models import MinIONEvent
from reads.models import MinIONEventType
from reads.models import MinIONRun
from reads.models import MinIONRunStats
from reads.models import MinIONRunStatus
from reads.models import MinIONScripts
from reads.models import MinIONStatus
from reads.models import MinIONmessages
from reads.models import RunSummary
from reads.models import JobMaster
from reads.models import ChannelSummary
from reads.models import HistogramSummary



class HistogramSummarySerializer(serializers.ModelSerializer):
    class Meta:
        model = HistogramSummary
        fields = (
            'id',
            'run_id',
            'read_type',
            'bin_width',
            'read_count',
            'read_length',
        )
        read_only = ('id',)


class ChannelSummarySerializer(serializers.ModelSerializer):
    class Meta:
        model = ChannelSummary
        fields = (
            'id',
            'run_id',
            'channel_number',
            'read_count',
            'read_length',
        )
        read_only = ('id',)


class MinIONRunStatsSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinIONRunStats
        fields = (
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
                  'start_time', 'run_id', 'type','created_date')


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


class RunSummarySerializer(serializers.HyperlinkedModelSerializer):
    typename = serializers.ReadOnlyField(source="type.name")

    class Meta:
        model = RunSummary
        fields = ('url', 'total_length', 'read_count', 'type', 'typename', 'max_length', 'min_length', 'pass_length', 'pass_max_length', 'pass_min_length', 'pass_count', 'id')
        read_only = ('id',)


class RunSummaryBarcodeSerializer(serializers.HyperlinkedModelSerializer):
    typename = serializers.ReadOnlyField(source="type.name")

    class Meta:
        model = RunSummaryBarCode
        fields = ('url', 'total_length', 'read_count', 'type', 'typename', 'max_length', 'min_length', 'pass_length', 'pass_max_length', 'pass_min_length', 'pass_count', 'id', 'barcode')
        read_only = ('id',)


class RunStatisticSerializer(serializers.ModelSerializer):
    typename = serializers.ReadOnlyField(source="type.name")

    class Meta:
        model = RunStatistic
        fields = ('total_length', 'read_count', 'type', 'typename', 'max_length', 'min_length', 'pass_length', 'pass_max_length', 'pass_min_length', 'pass_count', 'id', 'sample_time')
        read_only = ('id',)


class RunStatisticBarcodeSerializer(serializers.ModelSerializer):
    typename = serializers.ReadOnlyField(source="type.name")

    class Meta:
        model = RunStatisticBarcode
        fields = ('total_length', 'read_count', 'type', 'typename', 'max_length', 'min_length', 'pass_length', 'pass_max_length', 'pass_min_length', 'pass_count', 'id', 'sample_time', 'barcode')
        read_only = ('id',)

class JobSerializer(serializers.ModelSerializer):
    typename = serializers.ReadOnlyField(source="type.name")

    class Meta:
        model = JobMaster
        fields = ('run_id','job_name','var1','var2','var3','complete')



class MinIONRunSerializer(serializers.HyperlinkedModelSerializer):

    jobstodo = JobSerializer(
        many=True,
        read_only=True,
    )

    class Meta:
        model = MinIONRun
        fields = ('url', 'sample_name', 'minKNOW_version', 'minKNOW_flow_cell_id', 'run_name', 'run_id', 'is_barcoded','minION', 'barcodes', 'id', 'last_read', 'last_entry','jobstodo')
        read_only = ('id','sample_name','minKNOW_version', 'minKNOW_flow_cell_id', 'barcodes', 'last_read', 'last_entry','jobstodo')

