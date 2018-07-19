from datetime import datetime

from rest_framework import serializers

from devices.models import Flowcell, MinION
from jobs.models import JobMaster, JobType
from reads.models import (Barcode, BarcodeGroup, FastqRead, FastqReadExtra,
                          FastqReadType, FlowCellRun, MinIONControl,
                          MinIONEvent, MinIONEventType, MinionMessage,
                          MinIONRunStats, MinIONRunStatus, MinIONScripts,
                          MinIONStatus, Run, UserOptions, ChannelSummary, HistogramSummary,
                          RunStatisticBarcode, RunSummaryBarcode, GroupRun, FlowcellSummaryBarcode)


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
        'url', 'id', 'name', 'minION_name', 'space_available', 'minKNOW_version', 'status', 'computer', 'sample_name', 'run_status',
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


class MinionMessageSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinionMessage
        fields = (
        'minion', 'run', 'message', 'severity', 'timestamp', 'identifier',)
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

    sequence = serializers.CharField(allow_blank=True)
    quality = serializers.CharField(allow_blank=True)

    class Meta:
        model = FastqRead
        fields = (
            'url',
            'read_id',
            'read',
            'channel',
            'barcode',
            'sequence_length',
            'quality_average',
            'sequence',
            'quality',
            'is_pass',
            'start_time',
            'run',
            'type',
            'created_date'
        )

    def create(self, validated_data):

        #
        # The next two lines of code truncate the datetime information.
        # The RunStatisticsBarcode aggregates the date on minute level.
        #
        #start_time = validated_data['start_time']
        #start_time_truncated = datetime(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute)

        fastqread = FastqRead(
            read_id=validated_data['read_id'],
            read=validated_data['read'],
            channel=validated_data['channel'],
            barcode=validated_data['barcode'],
            sequence_length=validated_data['sequence_length'],
            quality_average=validated_data['quality_average'],
            is_pass=validated_data['is_pass'],
            start_time=validated_data['start_time'],
            run=validated_data['run'],
            type=validated_data['type']
        )

        fastqread.save()

        if fastqread.run.has_fastq:

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
        fields = ('id', 'minION', 'identifier', 'name', 'experiment_type', 'base_calling', 'flow_cell', 'kit', 'event_ratio', 'experiment_time')
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
            'experiment_type',
            'experiment_id',
            'fast5_output_fastq_in_hdf',
            'fast5_raw',
            'fast5_reads_per_folder',
            'fastq_enabled',
            'fastq_reads_per_file',
            'filename',
            'flowcell_type',
            'kit_classification',
            'local_basecalling',
            'sample_frequency',
            'sequencing_kit',
            'user_filename_input',

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
            'is_pass',
            'id',
            'barcode',
            'barcodename',
            'channel_presence'
            # 'number_active_channels'
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


class JobTypeSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:

        model = JobType
        fields = ('id', 'url', 'name', 'description', 'long_description', 'reference', 'transcriptome', 'readcount', 'private')


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

    minion_name = serializers.ReadOnlyField(source="minion.name")

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
            'minion_name',
            'id',
            'last_read',
            'last_entry',
            'max_channel',
            'flowcell_type',
            'jobstodo',
            'barcodes',
            'to_delete',
            'minion',
            'computer_name',
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
            'barcodes',
            'computer_name',
        )


class RunHistogramSummarySerializer(serializers.ModelSerializer):
    read_type_name = serializers.ReadOnlyField(source="read_type.name")
    barcode_name = serializers.ReadOnlyField(source="barcode.name")

    class Meta:
        model = HistogramSummary

        fields = (
            'barcode_name',
            'read_type_name',
            'bin_index',
            'read_count',
            'read_length'
        )


class FlowCellSerializer(serializers.HyperlinkedModelSerializer):

    barcodes = BarcodeSerializer(many=True, read_only=False, required=False)
    runs = RunSerializer(many=True, read_only=False, required=False)

    class Meta:
        model = Flowcell
        fields = ('url', 'name', 'runs', 'barcodes', 'id')

    def create(self, validated_data):

        # flowcell = Flowcell(
        #
        #     name=validated_data['name'],
        #     owner=validated_data['owner']
        # )
        # flowcell.save()

        flowcell = Flowcell(**validated_data)
        flowcell.save()

        job_type = JobType.objects.filter(name="ChanCalc")

        JobMaster.objects.create(
            flowcell=flowcell,
            job_type=job_type[0],
            last_read=0
        )

        return flowcell


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


class GroupRunSerializer(serializers.HyperlinkedModelSerializer):

    runs = RunSerializer(
        many=True,
        read_only=True
    )

    class Meta:

        model = GroupRun
        fields = ('id', 'url', 'name', 'device', 'runs')


class GroupRunMembershipSerializer(serializers.Serializer):

    grouprun_id = serializers.IntegerField()
    run_id = serializers.IntegerField()

    # grouprun = serializers.HyperlinkedRelatedField(
    #     view_name='grouprun-detail',
    #     queryset=GroupRun.objects.all()  # we should include a filter here
    # )

    # class Meta:
    #
    #     model = Run
    #     fields = ('url', 'grouprun')

    def create(self, validated_data):

        print('>>> validated_data')
        print(validated_data)
        print('<<< validated_data')
        #self.run.groupruns.add(self.grouprun)

        return GroupRun.objects.all()[0]



    # class Meta:
    #
    #     fields = ('grouprun', 'run')
    #


class FlowcellSummaryBarcodeSerializer(serializers.ModelSerializer):

    class Meta:
        model = FlowcellSummaryBarcode

        fields = (
            'flowcell',
            'read_type_name',
            'barcode_name',
            'status',
            'quality_sum',
            'read_count',
            'total_length',
            'max_length',
            'min_length',
            'channel_count',
            'average_read_length',
        )

        read_only = ('id',)

