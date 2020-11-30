import datetime
from pathlib import Path

from django.conf import settings
from rest_framework import serializers

from communication.models import NotificationConditions
from metagenomics.models import CentrifugeOutput
from minknow_data.models import Run, Flowcell
from reads.models import (Barcode, FastqFile, FastqRead,
                          FastqReadType, UserOptions, ChannelSummary, HistogramSummary,
                          RunStatisticBarcode, RunSummaryBarcode, FlowcellSummaryBarcode, JobType, JobMaster)


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


class FastqReadTypeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqReadType
        fields = ('url', 'name', 'id')
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
        read_only = (
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


class FastqReadNameSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqRead
        fields = ('read_id',)


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
        # Create tasks that need to be auto-created
        if created:
            JobMaster.objects.create(
                flowcell=flowcell,
                job_type=JobType.objects.get(name="UpdateFlowcellDetails"),
                last_read=0
            )
            # If user has allowed us to tweet them
            if UserOptions.objects.get(owner=flowcell.owner).tweet:
                # Auto create a Warnings and error notification
                NotificationConditions(notification_type="waer", flowcell=flowcell,
                                       creating_user=flowcell.owner).save()
                # Auto create a mux notification
                NotificationConditions(notification_type="mux", flowcell=flowcell,
                                       creating_user=flowcell.owner).save()
        return flowcell


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
            "read_count",
            "running",
            "complete",
            "paused",
            "server_initiated",
            "from_database"
        )


class JobMasterInsertSerializer(serializers.ModelSerializer):
    """
    Serialiser to validate and insert the JobMasters into the database, upon job creation via the reads/views.py
    get_or_create_jobs function
    """
    class Meta:
        """
        Meta options about model
        """
        model = JobMaster
        fields = (
            "job_type",
            "flowcell",
            "reference",
            "target_set",
            "from_database"
        )

    def validate(self, data):
        """
        Check that reference exists if type is Minimap2
        """
        if data['job_type'].id in [4, 16] and data['reference'] is None:
            raise serializers.ValidationError("Reference is mandatory for this task.")

        # If the job type is minimap check we can find the minimap2 executable
        if data['job_type'].id == 4:
            if not Path.exists(Path(getattr(settings, "MINIMAP2", None))):
                raise serializers.ValidationError("Minimap2 executable not found."
                                                  " Please configure in the environmental variables file.")

        # If it is not metagenomics double check the target set is None
        if data["job_type"].id != 10:
            data["target_set"] = None
        # Check if a metagenomics task is already running on this flowcell, if so refuse to create a new one
        # elif JobMaster.objects.filter(
        #         flowcell=data["flowcell"]).filter(job_type__name="Metagenomics").count() > 0:
        #     raise serializers.ValidationError("A metagenomics task is already running for this flowcell.")
        # If the job isn't EB or minimap2 or artic
        if data["job_type"].id not in [4, 15, 16]:
            data["reference"] = None
        if data["job_type"].id == 16:
            if JobMaster.objects.filter(
                flowcell=data["flowcell"]).filter(job_type_id=16).count() > 0:
                raise serializers.ValidationError("An artic task is already running for this flowcell.")
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
        :param validated_data: Sent validated form data from MinotourTask Controller to jobs.views task_list
        :return:
        """
        # Create a JobMaster
        job_master, created = JobMaster.objects.get_or_create(**validated_data)
        job_master.save()
        # If we have created a new job
        if created and not job_master.job_type_id == 17:
            # Get the flowcell that the JobMaster is for
            flowcell_to_update_activity = validated_data["flowcell"]
            # Update the last activity date
            flowcell_to_update_activity.last_activity_date = datetime.datetime.now(datetime.timezone.utc)
            # save the changes
            flowcell_to_update_activity.save()
        return job_master
