from django.conf import settings
from django.db import models


class ReferenceInfo(models.Model):

    reference_name = models.CharField(
        max_length=512
    )

    filename = models.CharField(
        max_length=512
    )

    bwa_index_file_location = models.CharField(
        max_length=512,
        null=True
    )

    minimap2_index_file_location = models.CharField(
        max_length=512,
        null=True
    )

    totalrefleN = models.BigIntegerField(

    )

    private = models.BooleanField(
        default=False
    )

    users = models.ManyToManyField(
        settings.AUTH_USER_MODEL,
        blank=True
    )

    def __str__(self):
        return self.reference_name


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
