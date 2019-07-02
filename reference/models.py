from django.conf import settings
from django.db import models


class ReferenceInfo(models.Model):
    """
    This class represents the genome references
    used in the mapping and metagenomics tasks.
    minoTour comes with the Ecoli reference.
    """

    name = models.CharField(

        max_length=512
    )

    filename = models.CharField(

        max_length=512
    )

    minimap2_index_file_location = models.CharField(

        max_length=512,
        null=True
    )

    length = models.BigIntegerField(

        default=0
    )

    private = models.BooleanField(

        default=False
    )

    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        related_name='user_references',
        on_delete=models.CASCADE,
        null=True
    )

    users = models.ManyToManyField(

        settings.AUTH_USER_MODEL,
        blank=True
    )

    transcripts = models.BooleanField(

        default=False
    )

    def __str__(self):
        return self.name


class ReferenceLine(models.Model):

    reference = models.ForeignKey(
        ReferenceInfo,
        on_delete=models.CASCADE,
        related_name='referencelines'
    )

    line_name = models.CharField(
        max_length=512
    )

    chromosome_length = models.BigIntegerField(

    )

    def __str__(self):
        return "{} {}".format(self.reference, self.line_name)
