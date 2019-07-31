"""
Django rest framework serialiser for the read until class
"""
from rest_framework.serializers import ListSerializer, Serializer
from readuntil.models import RejectedFastqRead, AcceptedFastqRead


class RejectedListSerializer(ListSerializer):
    """
    A seperate List serialiser class that will be invoked in the Rejected fastq read serialiser
    class when many is set to True, this class customises the behaviour of multiple creation to bulk insert
    """

    def create(self, validated_data):
        """
        Override the create function so that we bulk create
        :param validated_data: Validated data passed through a serialiser
        :return:
        """
        read_ids = [RejectedFastqRead(**item) for item in validated_data]
        return RejectedFastqRead.objects.bulk_create(read_ids)


class AcceptedListSerializer(ListSerializer):
    """
    A separate List serialiser class that will be invoked in the Accepted fastq read serialiser
    class when many is set to True, this class customises the behaviour of multiple creation to bulk insert
    """
    def create(self, validated_data):
        """
        Override the create function so that we bulk create
        :param validated_data: Validated data passed through a serialiser
        :return:
        """
        read_ids = [AcceptedFastqRead(**item) for item in validated_data]
        return AcceptedFastqRead.objects.bulk_create(read_ids)


class RejectedFastqSerialiser(Serializer):
    """
    A serialiser class for creating or deserialising the Rejected Fastq read objects in a read until experiment
    """

    class Meta:
        list_serializer_class = RejectedListSerializer


class AcceptedFastqSerialiser(Serializer):
    """
    A serialiser class for creating or deserialising the Accepted Fastq read objects in a read until experiment
    """

    class Meta:
        list_serializer_class = AcceptedListSerializer
