from django.db import models

from minknow_data.models import Run


# Create your models here.
class TomlFile(models.Model):
    """
    Store toml files as they start from an experiment, stored as JSON that needs to be translated back into TOML
    """
    run = models.ForeignKey(
        Run,
        related_name="readfish_run",
        on_delete=models.CASCADE
    )
    toml_json = models.TextField()
    sha_hash = models.CharField(
        max_length=64
    )


class CompletedBarcodes(models.Model):
    """
    Store Barcodes that have completed during the course of an Artic Run, so we skip them in our view code
    """
    barcode = models.CharField(
       max_length=10
    )
    run_id = models.CharField(
        max_length=64
    )

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['barcode', 'run_id'], name='unique completed barcode')
        ]