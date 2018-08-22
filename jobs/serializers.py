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
