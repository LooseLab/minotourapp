from rest_framework import serializers
from centrifuge.models import CentOutput, MappingTargets, MappingResults


class CentSerialiser(serializers.ModelSerializer):
    """
        :purpose: Serialise and deserialise Django Model objects from the CentOutput model.
        Used in centrifuge class
        :author: Rory
    """
    class Meta:
        model = CentOutput
        fields = ("name", "tax_id")


class MappedResultsSerialiser(serializers.ModelSerializer):
    """Serialiser to represent CartographyMapped model TODO currently unused"""
    class Meta:
        model = MappingResults
        fields = ("species", "tax_id", "total_reads", "alert_level", "num_reads", "red_reads", "task_meta")


class MappedTargetsSerialiser(serializers.ModelSerializer):
    """Serialiser to represent the CartographyGuide model TODO currently unsused"""
    class Meta:
        model = MappingTargets
        fields = ("species", "tax_id", "cart_mapped_id")


