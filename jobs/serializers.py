from rest_framework import serializers

from jobs.models import JobMaster, JobType
from centrifuge.models import MappingTarget
from reference.models import ReferenceInfo


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

    def create(self, validated_data):
        """
        Create the mapping task if the reference type is metagenoimics
        :param validated_data:
        :return:
        """

        job_master, created = JobMaster.objects.get_or_create(**validated_data)
        job_master.save()
        checklist = []
        if created and validated_data["job_type"].name == "Metagenomics":

            job_type = JobType.objects.get(name="Other")
            reference = list(MappingTarget.objects.all().values_list("species", "name", "gff_line_type"))
            reference.append(("Bacillus_anthracis", "", "gene"))
            for ref in reference:
                init = ref[0] + ref[2]

                if init in checklist and ref[2] != "plasmid":
                    continue
                else:
                    checklist.append(init)

                if ref[2] == "plasmid":
                    refer = ReferenceInfo.objects.get(name=ref[0].replace(" ", "_") + "_" + ref[1])

                else:
                    refer = ReferenceInfo.objects.get(name=ref[0].replace(" ", "_"))

                jm = JobMaster(job_type=job_type, reference=refer, flowcell=validated_data["flowcell"],
                               complete=True)
                jm.save()

        return job_master

