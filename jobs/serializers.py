from rest_framework import serializers

from jobs.models import JobMaster, JobType
from centrifuge.models import MappingTarget, CentrifugeOutput
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
        # Check if a metagenomics task is already running on this flowcell, if so refuse to create a new one
        if data["job_type"].id == 10 and JobMaster.objects.filter(
                flowcell=data["flowcell"]).filter(job_type__name="Metagenomics").count() > 0:

            raise serializers.ValidationError("A metagenomics task is already running for this flowcell.")
        # If there's a target set
        if data["job_type"].id == 10:
            data["reference"] = None
            data["target_set"] = MappingTarget.objects.values_list("target_set", flat=True).distinct()[0]

        if data["job_type"].id == 13:
            metagenomics_task = JobMaster.objects.filter(
                flowcell=data["flowcell"]).filter(job_type__name="Metagenomics")
            if metagenomics_task.count() != 1:
                raise serializers.ValidationError("No metagenomics task found for this flowcell,"
                                                  " please start a metagenomics task before calculating"
                                                  " a lineages diagram")
            elif metagenomics_task.count() > 0 and CentrifugeOutput.objects.filter(task=metagenomics_task).count() < 10:
                raise serializers.ValidationError("Not enough metagenomics data, please try again in a few minutes")

        return data

    def create(self, validated_data):
        """
        Create a new job master, courtesy of jobs./views.py task_list post function
        :param validated_data: Sent validated form data from MTask Controller to jobs.views task_list
        :return:
        """

        # Create a JobMaster
        job_master, created = JobMaster.objects.get_or_create(**validated_data)
        job_master.save()
        # If we have created a new job and it's a metagenomics job
        if created and validated_data["job_type"].name == "Metagenomics":
            # Set the new target set
            # Get the jobtype for mapping target reads
            job_type = JobType.objects.get(name="Other")
            # get the target set we need
            mapping_target = list(MappingTarget.objects.filter(target_set=validated_data["target_set"])
                                  .values_list("species", "name", "gff_line_type"))
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
                    refer = ReferenceInfo.objects.get(name=hold[0].replace(" ", "_") + "_" + hold[1])
                else:
                    refer = ReferenceInfo.objects.get(name=hold[0].replace(" ", "_"))

                jm, created = JobMaster.objects.get_or_create(job_type=job_type, reference=refer,
                                                     flowcell=validated_data["flowcell"],
                                                     complete=True)
                jm.save()

        return job_master
