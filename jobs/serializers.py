from rest_framework import serializers

from jobs.models import JobMaster


class JobMasterSerializer(serializers.ModelSerializer):

    task_type_name = serializers.ReadOnlyField(

        source="job_type.name"
    )

    reference_name = serializers.ReadOnlyField(

        source="reference.name"
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
            "complete"
        )


class JobMasterInsertSerializer(serializers.ModelSerializer):

    class Meta:

        model = JobMaster

        fields = (
            "job_type",
            "flowcell",
            "reference"
        )

    def validate(self, data):
        """
        Check that reference exists if type is Minimap2
        """

        if data['job_type'].id == 4 and data['reference'] is None:

            raise serializers.ValidationError("Reference is mandatory for task Minimap2")
        # CHeck if a metagenomics task is already running on this flowcell, if so refuse to create a new one
        if data["job_type"].id == 10 and JobMaster.objects.filter(
                flowcell=data["flowcell"]).filter(job_type__name="Metagenomics").count() > 0:

            raise serializers.ValidationError("A metagnomics task is already running for this flowcell.")

        return data
