from rest_framework import serializers
from centrifuge.models import CentOutput, CartographyGuide, CartographyMapped, LineageValues

class CentSerialiser(serializers.ModelSerializer):
    """
        :purpose: Serialise and deserialise Django Model objects from the CentOutput model.
        Used in centrifuge class
        :author: Rory
    """
    class Meta:
        model = CentOutput
        fields = ("name", "tax_id", "num_matches", "sum_unique", "taxa")


class CartMappedSerialiser(serializers.ModelSerializer):
    """Serialiser to represent CartographyMapped model TODO currently unused"""
    class Meta:
        model = CartographyMapped
        fields = ("species", "tax_id", "total_reads", "alert_level", "num_reads", "red_reads", "task_meta")


class CartGuideSerialiser(serializers.ModelSerializer):
    """Serialiser to represent the CartographyGuide model TODO currently unsused"""
    class Meta:
        model = CartographyGuide
        fields = ("species", "tax_id", "cart_mapped_id")

