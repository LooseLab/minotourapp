"""Model classes for the artic app."""
from django.db import models

from minknow_data.models import Flowcell
from reads.models import JobMaster, Barcode


class ArticFireConditions(models.Model):
    """
    Shows the conditions for each task that wil cause the Fire Artic JobMaster to be created
    """
    flowcell = models.ForeignKey(
        Flowcell,
        related_name="flowcell_firing_conditions",
        on_delete=models.CASCADE,
    )
    percent_of_amplicons = models.IntegerField(default=90)
    x_coverage = models.IntegerField(default=20)
    percent_amp_at_x_cov = models.FloatField(default=0.0)


class ArticBarcodeMetadata(models.Model):
    """
    Contains the metadata for each barcode in an artic track coverage task.
    """
    average_coverage = models.IntegerField(default=0)
    variance_coverage = models.FloatField(default=0.0)
    minimum_coverage = models.IntegerField(default=0)
    maximum_coverage = models.IntegerField(default=0)
    percentage_of_reads_in_barcode = models.FloatField(default=0)
    has_sufficient_coverage = models.BooleanField(default=False)
    has_finished = models.BooleanField(default=False)
    marked_for_rerun = models.BooleanField(default=False)
    projected_to_finish = models.BooleanField(default=False)
    final_completion = models.BooleanField(default=False,null=True)
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
    artic_fire_conditions = models.ManyToManyField(ArticFireConditions)

    def __str__(self):
        return f"{self.flowcell} {self.job_master} {self.barcode}"
