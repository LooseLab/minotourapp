import datetime
from pathlib import Path

from django.conf import settings
from rest_framework import serializers

from centrifuge.models import CentrifugeOutput, MappingTarget
from reads.models import (Barcode, FastqFile, FastqRead,
                          FastqReadType, MinIONControl,
                          MinIONEvent, MinIONEventType, MinionMessage,
                          MinIONRunStats, MinIONRunStatus, MinIONScripts,
                          MinIONStatus, Run, UserOptions, ChannelSummary, HistogramSummary,
                          RunStatisticBarcode, RunSummaryBarcode, GroupRun, FlowcellSummaryBarcode, Flowcell, MinION,
                          FlowcellTab, JobType, JobMaster)
from reference.models import ReferenceInfo


class FastqFileSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqFile
        fields = (
            'url',
            'id',
            'name',
            'runid',
            'md5',
        )
        read_only=('id',)


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


#class MinIONRunStatsSerializer(serializers.HyperlinkedModelSerializer):
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
            'pore',
            'no_pore',
            'zero',
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
        fields = ('url', 'name', 'id')
        read_only = ('id',)


class MinionMessageSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinionMessage
        fields = (
        'minion', 'run', 'message', 'severity', 'timestamp', 'identifier', 'full_text')
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


class FastqReadGetSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = FastqRead
        fields = (
            'url',
            'read_id',
            'read',
            'channel',
            'barcode',
            'barcode_name',
            'sequence_length',
            'quality_average',
            'is_pass',
            'start_time',
            'run',
            'type',
            'created_date'
        )


class FastqReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = FastqRead
        fields = (
            'read_id',
            'read',
            'channel',
            'barcode',
            'barcode_name',
            'rejected_barcode',
            'sequence_length',
            'quality_average',
            'sequence',
            'quality',
            'is_pass',
            'start_time',
            'run',
            'type',
            'created_date',
            'fastqfile'
        )

    def create(self, validated_data):

        run = validated_data['run']

        fastqread = FastqRead(
            read_id=validated_data['read_id'],
            read=validated_data['read'],
            channel=validated_data['channel'],
            barcode=validated_data['barcode'],
            rejected_barcode=validated_data['rejected_barcode'],
            barcode_name=validated_data['barcode_name'],
            sequence_length=validated_data['sequence_length'],
            quality_average=validated_data['quality_average'],
            sequence=validated_data['sequence'],
            quality=validated_data['quality'],
            is_pass=validated_data['is_pass'],
            start_time=validated_data['start_time'],
            run=run,
            flowcell=run.flowcell,
            type=validated_data['type'],
            fastqfile=validated_data['fastqfile']
        )

        fastqread.save()

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

    class Meta:
        model = Barcode

        fields = (
            'url',
            'id',
            'name',
            'run',
        )

        read_only = (
            'url',
            'id'
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


class FlowcellSerializer(serializers.HyperlinkedModelSerializer):

    barcodes = BarcodeSerializer(many=True, read_only=False, required=False)

    runs = RunSerializer(many=True, read_only=False, required=False)



    class Meta:

        model = Flowcell

        fields = ('url', 'name', 'runs', 'barcodes', 'id', 'has_fastq')

    def create(self, validated_data):

        #
        # Flowcell names must be unique for a particular user
        #
        flowcell, created = Flowcell.objects.get_or_create(**validated_data)
        flowcell.save()

        if created:

            job_type = JobType.objects.filter(name="ChanCalc")

            JobMaster.objects.create(
                flowcell=flowcell,
                job_type=job_type[0],
                last_read=0
            )

            JobMaster.objects.create(
                flowcell=flowcell,
                job_type=JobType.objects.filter(name="UpdateFlowcellDetails")[0],
                last_read=0
            )

        return flowcell


class GroupRunSerializer(serializers.HyperlinkedModelSerializer):

    runs = RunSerializer(
        many=True,
        read_only=True
    )

    class Meta:

        model = GroupRun
        fields = ('id', 'url', 'name', 'device', 'runs')


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


class JobMasterSerializer(serializers.ModelSerializer):
    """
    Serialise the JobMaster django model, used to fetch JobMasters
    """
    task_type_name = serializers.ReadOnlyField(
        source="job_type.name"
    )

    reference_name = serializers.ReadOnlyField(
        default="N/A",
        source="reference.name"
    )
    server_initiated = serializers.ReadOnlyField(
        source="job_type.private"
    )

    class Meta:

        model = JobMaster

        fields = (
            "id",
            "flowcell",
            "task_type_name",
            "reference_name",
            "last_read",
            "read_count",
            "running",
            "complete",
            "paused",
            "server_initiated"
        )


class JobMasterInsertSerializer(serializers.ModelSerializer):
    """
    Serialiser to validate and insert the JobMasters into the database, upon job creation via the reads/views.py
    get_or_create_jobs function
    """
    class Meta:

        model = JobMaster

        fields = (
            "job_type",
            "flowcell",
            "reference",
            "target_set"
        )

    def validate(self, data):
        """
        Check that reference exists if type is Minimap2
        """
        if data['job_type'].id == 4 and data['reference'] is None:

            raise serializers.ValidationError("Reference is mandatory for task Minimap2")

        # If the job type is minimap check we can find the minimap2 executable
        if data['job_type'].id == 4:
            if not Path.exists(Path(getattr(settings, "MINIMAP2", None))):
                raise serializers.ValidationError("Minimap2 executable not found."
                                                  " Please configure in the environmental variables file.")

        # If it is not metagenomics double check the target set is None
        if data["job_type"].id != 10:

            data["target_set"] = None

        # Check if a metagenomics task is already running on this flowcell, if so refuse to create a new one
        elif JobMaster.objects.filter(
                flowcell=data["flowcell"]).filter(job_type__name="Metagenomics").count() > 0:

            raise serializers.ValidationError("A metagenomics task is already running for this flowcell.")

        # If the job isn't EB or minimap2
        if data["job_type"].id not in [4, 15]:

            data["reference"] = None

        # If the job is to calculate a sankey diagram
        if data["job_type"].id == 13:
            # Get the metagenomics task object
            metagenomics_task = JobMaster.objects.filter(
                flowcell=data["flowcell"]).filter(job_type__name="Metagenomics")

            # If we don't have a metagenomics task we can't calculate a sankey diagram
            if metagenomics_task.count() != 1:

                raise serializers.ValidationError("No metagenomics task found for this flowcell,"
                                                  " please start a metagenomics task before calculating"
                                                  " a lineages diagram")

            # If we don't have enough data to generate a sankey
            elif CentrifugeOutput.objects.filter(task=metagenomics_task[0]).count() < 10:

                raise serializers.ValidationError("Not enough metagenomics data, please try again in a few minutes")

        return data

    def create(self, validated_data):
        """
        Create a new job master, courtesy of jobs./views.py ge_or_create_tasks post function
        :param validated_data: Sent validated form data from MTask Controller to jobs.views task_list
        :return:
        """
        print(validated_data)
        # Create a JobMaster
        job_master, created = JobMaster.objects.get_or_create(**validated_data)
        job_master.save()
        print(created)
        # If we have created a new job
        if created:
            # Get the flowcell that the JobMaster is for
            flowcell_to_update_activity = validated_data["flowcell"]
            print(flowcell_to_update_activity)
            # Update the last activity date
            flowcell_to_update_activity.last_activity_date = datetime.datetime.now(datetime.timezone.utc)
            # save the changes
            flowcell_to_update_activity.save()

            # If we have created a metagenomics task
            if validated_data["job_type"].name == "Metagenomics":
                # We have to create "Other" job types to do the validation mapping
                # Get the jobtype for mapping target reads
                job_type = JobType.objects.get(name="Other")

                # get the target set we need
                if "target_set" is not None:
                    mapping_target = list(MappingTarget.objects.filter(target_set=validated_data["target_set"])
                                          .values_list("species", "name", "gff_line_type"))

                else:
                    mapping_target = None

                if mapping_target:
                    # Only create one gene entry per species, by counting how many are created in this dictionary,
                    # doesn't really need to be a dict
                    counts_dict = {}
                    # Set to store unique species, target type ( gene or plasmid )
                    setty = set()
                    # Holding copy is a copy of the mapping targets, to append to below
                    holding_copy = mapping_target
                    # loop through references and create a job to map any that are only plasmid targets have
                    for target in mapping_target:
                        # Get unique types of gff regions for each species
                        setty.add((target[0], target[2]))
                        # If we don't have just for the reference but have one for a plasmid,
                        # add a new gene entry to holding copy
                        if (target[0], "gene") not in setty and (target[0], "", "gene") not in mapping_target:
                            holding_copy.append((target[0], "", "gene"))
                    # Loop through and create the mapping job masters
                    for hold in holding_copy:
                        # Check to see whether we have a gene entry for this species in the dictionary, if we do skip
                        if (hold[0], hold[2]) not in counts_dict and hold[2] == "gene":
                            counts_dict[(hold[0], hold[2])] = 1
                        elif (hold[0], hold[2]) in counts_dict and hold[2] == "gene":
                            continue
                        # get the correct reference for that JobMaster
                        if hold[2] == "plasmid":
                            # TODO temporary fix
                            refer = ReferenceInfo.objects.filter(name=hold[0].replace(" ", "_") + "_" + hold[1]).first()
                        else:
                            # TODO temporary fix
                            refer = ReferenceInfo.objects.filter(name=hold[0].replace(" ", "_")).first()

                        jm, created = JobMaster.objects.get_or_create(job_type=job_type, reference=refer,
                                                             flowcell=validated_data["flowcell"],
                                                             complete=True)
                        jm.save()
            # ############################ End of other mapping creation
        return job_master