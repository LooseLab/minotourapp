"""Model classes for the artic app."""
from django.db import models

from reads.models import Flowcell, JobMaster, Barcode


class ArticBarcodeMetadata(models.Model):
    """
    Contains the metadata for each barcode in an artic track coverage task.
    """
    average_coverage = models.IntegerField(default=0)
    variance_coverage = models.FloatField(default=0.0)
    minimum_coverage = models.IntegerField(default=0)
    maximum_coverage = models.IntegerField(default=0)
    percentage_of_reads_in_barcode = models.FloatField(default=0)
    percentage_bases_over_20x = models.FloatField(default=0.0)
    percentage_bases_over_250x = models.FloatField(default=0.0)
    percentage_bases_over_200x = models.FloatField(default=0.0)
    has_sufficient_coverage = models.BooleanField(default=False)
    has_finished = models.BooleanField(default=False)
    marked_for_rerun = models.BooleanField(default=False)
    flowcell = models.ForeignKey(
        Flowcell,
        related_name="FlowcellsArticBarcodeMetadatas",
        on_delete=models.CASCADE,
        null=True
    )
    job_master = models.ForeignKey(
        JobMaster,
        related_name="JobMastersArticBarcodeMetadatas",
        on_delete=models.CASCADE,
        null=True
    )
    barcode = models.ForeignKey(
        Barcode,
        related_name="barcodeArticBarcodeMetadata",
        on_delete=models.CASCADE,
        null=True
    )

    def __str__(self):
        return f"{self.flowcell} {self.job_master} {self.barcode}"