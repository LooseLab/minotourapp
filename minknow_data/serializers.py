from rest_framework import serializers

from minknow_data.models import MinionMessage, MinionControl, MinionEventType, MinionEvent, MinionScripts, MinionInfo, \
    MinionRunInfo, MinionRunStats, Minion


class MinionMessageSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinionMessage
        fields = (
        'minion', 'run', 'message', 'severity', 'timestamp', 'identifier', 'full_text')
        read_only = ('id',)


class MinIONControlSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinionControl
        fields = ('id', 'minION', 'job', 'custom', 'complete')
        read_only = ('id',)


class MinIONEventTypeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinionEventType
        fields = ('url', 'name',)
        read_only = ('id',)


class MinIONEventSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinionEvent
        fields = ('id', 'computer_name', 'minION', 'event', 'datetime')
        read_only = ('id',)


class MinIONScriptsSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinionScripts
        fields = ('id', 'minION', 'identifier', 'name', 'experiment_type', 'base_calling', 'flow_cell', 'kit', 'event_ratio', 'experiment_time')
        read_only = ('id',)


class MinionInfoSerializer(serializers.HyperlinkedModelSerializer):
    """
    Serialiser for minion_info model.
    """
    class Meta:
        model = MinionInfo
        fields = (
            'id',
            'minion',
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
            'minknow_version'
        )
        read_only = ('id',)


class MinionRunInfoSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinionRunInfo
        fields = (
            'id',
            'minion',
            'minION_name',
            # 'minKNOW_status',
            'minKNOW_current_script',
            'minKNOW_sample_name',
            'minKNOW_exp_script_purpose',
            'minKNOW_flow_cell_id',
            'minKNOW_run_name',
            'minKNOW_version',
            'run',
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
            'wells_per_channel',
            'target_temp'
        )
        read_only = ('id', 'minION_name')


class MinionRunStatsSerializer(serializers.HyperlinkedModelSerializer):
    """
    Serialiser for MinION run stats. Stats about the Run provided by minFQ/minKNOW.
    """
    class Meta:
        model = MinionRunStats
        fields = (
            'id',
            'minion',
            'run',
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
            'pore',
            'no_pore',
            'zero',
            'occupancy',
            'minKNOW_read_count',
            'minKNOW_histogram_values',
            'minKNOW_histogram_bin_width',
            'actual_max_val',
            'n50_data',
            'estimated_selected_bases',
            'basecalled_bases',
            'basecalled_fail_read_count',
            'basecalled_pass_read_count'
        )
        read_only = ('id', 'occupancy',)


class MinionSerializer(serializers.HyperlinkedModelSerializer):
    """
    Serialiser for the Minion model.
    """
    channel_count = serializers.SerializerMethodField()

    class Meta:
        model = Minion
        fields = (
        'url', 'id', 'name', 'minION_name', 'space_available', 'status', 'computer', 'sample_name', 'run_status',
        'flow_cell_id', 'run_name', 'total_drive_space', 'space_till_shutdown', 'warnings', 'last_run', 'current_script',
        'event_yield', 'voltage_value', 'channel_count', 'flow_cell_id')
        read_only = (
        'status', 'space_available', 'computer', 'sample_name', 'run_status', 'flow_cell_id',
        'run_name', 'total_drive_space', 'space_till_shutdown', 'warnings', 'last_run', 'current_script', 'event_yield',
        'voltage_value',)

    def get_channel_count(self, obj):
        if hasattr(obj.currentrundetails.last(), "run_id"):
            return obj.currentrundetails.last().run.flowcell.size
        return "Unknown"