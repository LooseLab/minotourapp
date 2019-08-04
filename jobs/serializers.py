"""
Serialisers for the JobMaster and JobTypes django models
"""
from rest_framework import serializers
from jobs.models import JobMaster, JobType
from centrifuge.models import MappingTarget, CentrifugeOutput
from reference.models import ReferenceInfo
import datetime
from pathlib import Path
from django.conf import settings


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
            elif metagenomics_task.count() > 0 and CentrifugeOutput.objects.filter(
                    task=metagenomics_task)[:20].count() < 10:

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
                            refer = ReferenceInfo.objects.get(name=hold[0].replace(" ", "_") + "_" + hold[1])
                        else:
                            refer = ReferenceInfo.objects.get(name=hold[0].replace(" ", "_"))

                        jm, created = JobMaster.objects.get_or_create(job_type=job_type, reference=refer,
                                                             flowcell=validated_data["flowcell"],
                                                             complete=True)
                        jm.save()
            # ############################ End of other mapping creation
        return job_master
