from datetime import datetime

import numpy as np
from rest_framework import serializers

from devices.models import Flowcell, MinION
from jobs.models import JobMaster, JobType
from reads.models import (Barcode, BarcodeGroup, FastqRead, FastqReadExtra,
                          FastqReadType, FlowCellRun, MinIONControl,
                          MinIONEvent, MinIONEventType, MinIONmessages,
                          MinIONRunStats, MinIONRunStatus, MinIONScripts,
                          MinIONStatus, Run, UserOptions)
from stats.models import (ChannelSummary, HistogramSummary,
                          RunStatisticBarcode, RunSummaryBarcode)


class UserOptionsSerializer(serializers.ModelSerializer):
    class Meta:
        model = UserOptions
        fields = (
            'id',
            'owner',
            'twitterhandle',
            'tweet',
            'email',
        )
        read_only = ('id',)


class ChannelSummarySerializer(serializers.ModelSerializer):

    class Meta:

        model = ChannelSummary

        fields = (
            'channel',
            'read_count',
            'read_length',
        )


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
        read_only = ('id', 'occupancy',)


class MinIONSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinION
        fields = (
        'url', 'space_available', 'minKNOW_version', 'minION_name', 'status', 'computer', 'sample_name', 'run_status',
        'flow_cell_id', 'run_name', 'total_drive_space', 'space_till_shutdown', 'warnings', 'last_run', 'currentscript',
        'event_yield', 'voltage_value',)
        read_only = (
        'status', 'space_available', 'minKNOW_version', 'computer', 'sample_name', 'run_status', 'flow_cell_id',
        'run_name', 'total_drive_space', 'space_till_shutdown', 'warnings', 'last_run', 'currentscript', 'event_yield',
        'voltage_value',)


class FastqReadTypeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqReadType
        fields = ('url', 'name',)
        read_only = ('id',)


class MinIONmessagesSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinIONmessages
        fields = (
        'minION', 'run_id', 'minKNOW_message', 'minKNOW_severity', 'minKNOW_message_timestamp', 'minKNOW_identifier',)
        read_only = ('id',)


class MinIONControlSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinIONControl
        fields = ('id', 'minION', 'job', 'custom', 'complete')
        read_only = ('id',)


class MinIONEventTypeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinIONEventType
        fields = ('url', 'name',)
        read_only = ('id',)


class FastqReadSerializer(serializers.HyperlinkedModelSerializer):
    # class FastqReadSerializer(serializers.Serializer):

    sequence = serializers.CharField(allow_blank=True)
    quality = serializers.CharField(allow_blank=True)

    # sequence = serializers.ReadOnlyField(source="extra.sequence")
    # quality = serializers.ReadOnlyField(source="extra.quality")

    class Meta:
        model = FastqRead
        fields = (
        'url', 'read_id', 'read', 'channel', 'barcode', 'sequence_length', 'quality_average', 'sequence', 'quality',
        'is_pass',
        'start_time', 'run_id', 'type', 'created_date')

    def create(self, validated_data):

        #
        # The next two lines of code truncate the datetime information.
        # The RunStatisticsBarcode aggregates the date on minute level.
        #
        start_time = validated_data['start_time']
        start_time_truncated = datetime(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute)

        fastqread = FastqRead(
            read_id=validated_data['read_id'],
            read=validated_data['read'],
            channel=validated_data['channel'],
            barcode=validated_data['barcode'],
            # sequence_length = len(validated_data['sequence']),
            sequence_length=validated_data['sequence_length'],
            # quality_average = qualmean,
            quality_average=validated_data['quality_average'],
            is_pass=validated_data['is_pass'],
            start_time=start_time_truncated,
            run_id=validated_data['run_id'],
            type=validated_data['type']
        )

        # print('>>>> fastqread')
        # print(fastqread)

        # print('>>>> save')
        fastqread.save()
        if fastqread.run_id.has_fastq:
            fastqread_extra = FastqReadExtra(
                fastqread=fastqread,
                sequence=validated_data['sequence'],
                quality=validated_data['quality']
            )

            fastqread_extra.save()

        return fastqread


class FastqReadNameSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqRead
        fields = ('read_id',)


class MinIONEventSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinIONEvent
        fields = ('id', 'computer_name', 'minION', 'event', 'datetime')
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
            # 'minKNOW_status',
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
            'minKNOW_colours_string',
            'minKNOW_computer',
        )
        read_only = ('id', 'minION_name')


class RunSummaryBarcodeSerializer(serializers.HyperlinkedModelSerializer):
    typename = serializers.ReadOnlyField(
        source="type.name"
    )

    barcodename = serializers.ReadOnlyField(
        source="barcode.name"
    )

    class Meta:
        model = RunSummaryBarcode

        fields = (
            'url',
            'total_length',
            'read_count',
            'type',
            'typename',
            'max_length',
            'min_length',
            'pass_length',
            'pass_max_length',
            'pass_min_length',
            'pass_count',
            'id',
            'barcode',
            'barcodename',
            'number_active_channels'
        )

        read_only = ('id',)


class RunStatisticBarcodeSerializer(serializers.ModelSerializer):
    typename = serializers.ReadOnlyField(
        source="type.name"
    )

    barcodename = serializers.ReadOnlyField(
        source="barcode.name"
    )

    class Meta:
        model = RunStatisticBarcode

        fields = (
            'total_length',
            'read_count',
            'type',
            'typename',
            'quality_sum',
            'pass_quality_sum',
            'max_length',
            'min_length',
            'pass_length',
            'pass_max_length',
            'pass_min_length',
            'pass_count',
            'id',
            'sample_time',
            'barcode',
            'barcodename',
            'number_active_channels'
        )

        read_only = ('id',)


class JobTypeSerializer(serializers.ModelSerializer):
    # type_name = GroupRunTypeSerializer(read_only=True, many=True)

    class Meta:
        model = JobType
        fields = ('name', 'description', 'reference', 'readcount', 'transcriptome', 'type_name',)


class JobSerializer(serializers.ModelSerializer):
    typename = serializers.ReadOnlyField(source="type.name")

    class Meta:
        model = JobMaster
        fields = ('run_id', 'job_type', 'reference', 'last_read', 'tempfile_name', 'read_count', 'complete', 'running')


class BarcodeSerializer(serializers.HyperlinkedModelSerializer):
    barcodegroupname = serializers.ReadOnlyField(
        source="barcodegroup.name"
    )

    class Meta:
        model = Barcode

        fields = (
            'url',
            'id',
            'name',
            'barcodegroupname',
            'run',
        )

        read_only = (
            'url',
            'id'
        )


class BarcodeGroupSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = BarcodeGroup

        fields = (
            'id',
            'name',
            'flowcell',
            'url'
        )

        read_only = (
            'id',
            'url'
        )

class RunSerializer(serializers.HyperlinkedModelSerializer):

    jobstodo = JobSerializer(
        many=True,
        read_only=True,
    )

    barcodes = BarcodeSerializer(
        many=True,
        read_only=True
    )

    class Meta:
        model = Run
        fields = (
            'url',
            'flowcell_name',
            'flowcell',
            'start_time',
            'active',
            'sample_name',
            'minKNOW_version',
            'minKNOW_flow_cell_id',
            'name',
            'runid',
            'is_barcoded',
            'has_fastq',
            'minion',
            'id',
            'last_read',
            'last_entry',
            'max_channel',
            'flowcell_type',
            'jobstodo',
            'barcodes',
            'to_delete',
        )

        read_only = (
            'id',
            'flowcell_name',
            'flowcell',
            'start_time',
            'sample_name',
            'minKNOW_version',
            'minKNOW_flow_cell_id',
            'last_read',
            'last_entry',
            'jobstodo',
            'barcodes'
        )


class RunHistogramSummarySerializer(serializers.ModelSerializer):
    read_type_name = serializers.ReadOnlyField(source="read_type.name")
    barcode_name = serializers.ReadOnlyField(source="barcode.name")

    class Meta:
        model = HistogramSummary

        fields = (
            'barcode_name',
            'read_type_name',
            'bin_width',
            'read_count',
            'read_length'
        )


class FlowCellSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Flowcell
        fields = ('url', 'id', 'name',)


class FlowCellRunSerializer(serializers.HyperlinkedModelSerializer):
    barcodes = BarcodeSerializer(
        many=True,
        read_only=True
    )

    class Meta:
        model = FlowCellRun
        fields = (
        'flowcell', 'id', 'run_name', 'run', 'name', 'barcodes', 'last_entry', 'start_time', 'last_read', 'sample_name',
        'minKNOW_flow_cell_id', 'minKNOW_version', 'active')
        read_only = (
            'id',
            'run_name',
            'name',
            'barcodes',
            'last_entry', 'start_time', 'last_read', 'sample_name', 'minKNOW_flow_cell_id', 'minKNOW_version', 'active'
        )
