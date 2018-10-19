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

        print(data)

        print(data['job_type'].id)
        print(data['reference'])

        if data['job_type'].id == 4 and data['reference'] is None:

            raise serializers.ValidationError("Reference is mandatory for task Minimap2")

        return data
