from rest_framework import serializers
from centrifuge.models import CentOutput, CartographyGuide, CartographyMapped, LineageValues

class CentSerialiser(serializers.ModelSerializer):
    """Serialiser to represent CentOutput model"""
    class Meta:
        model = CentOutput
        fields = ("name", "tax_id", "num_matches", "sum_unique", "taxa")


class CartMappedSerialiser(serializers.ModelSerializer):
    """Serialiser to represent CartographyMapped model"""
    class Meta:
        model = CartographyMapped
        fields = ("species", "tax_id", "total_reads", "alert_level", "num_reads", "red_reads", "task_meta")


class CartGuideSerialiser(serializers.ModelSerializer):
    """Serialiser to represent the CartographyGuide model"""
    class Meta:
        model = CartographyGuide
        fields = ("species", "tax_id", "cart_mapped_id")

