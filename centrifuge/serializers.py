from rest_framework import serializers
from centrifuge.models import CentrifugeOutput


class CentListSerialiser(serializers.ModelSerializer):
    """
        :purpose: Serialise and deserialise Django Model objects from the CentOutput model.
        Used in centrifuge class
        :author: Rory
    """
    class Meta:
        model = CentrifugeOutput
        fields = '__all__'

        def create(self, validated_data):
            CentrifugeOutput(
                name=validated_data["name"],

                tax_id = validated_data["tax_id"],

                task = models.ForeignKey(
                    JobMaster,
                    related_name="centrifuge_output",
                    on_delete=models.CASCADE,
                ),

                num_matches = models.IntegerField(default=0),

                sum_unique = models.IntegerField(default=0),

                barcode_name = models.CharField(max_length=50),

                proportion_of_classified = models.FloatField(default=0.0),

                superkingdom = models.CharField(null=True, max_length=70),

                phylum = models.CharField(null=True, max_length=70),

                classy = models.CharField(null=True, max_length=70),

                order = models.CharField(null=True, max_length=70),

                family = models.CharField(null=True, max_length=70),

                genus = models.CharField(null=True, max_length=70),

                species = models.CharField(null=True, max_length=70)
            )


