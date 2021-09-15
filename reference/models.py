import os

from django.conf import settings
from django.db import models
from django.db.models.signals import post_delete
from django.dispatch import receiver


class ReferenceInfo(models.Model):
    """
    This class represents the genome references
    used in the mapping and metagenomics tasks.
    minoTour comes with the Ecoli reference.
    """

    name = models.CharField(
        max_length=256
    )
    file_name = models.CharField(
        max_length=256
    )
    minimap2_index_file_location = models.CharField(
        default="",
        max_length=256
    )
    file_location = models.FileField(
        upload_to="reference_files",
        default="",
        null=True,
        max_length=256
    )
    length = models.BigIntegerField(
        default=0
    )
    private = models.BooleanField(
        default=False
    )
    uploader = models.ForeignKey(
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
    sha256_checksum = models.CharField(
        max_length=256,
        default=""
    )

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=["name", "private", "uploader"], name='unique_reference'),
        ]

    def __str__(self):
        return f"Name: {self.name}, Private: {self.private}"



class ReferenceLine(models.Model):
    """
    Represents contigs inside the reference
    """
    reference = models.ForeignKey(
        ReferenceInfo,
        on_delete=models.CASCADE,
        related_name='reference_lines'
    )
    line_name = models.CharField(
        max_length=512
    )
    chromosome_length = models.BigIntegerField(
    )
    def __str__(self):
        return f"{self.reference} {self.line_name}"


@receiver(post_delete, sender=ReferenceInfo)
def delete_reference_files(sender, instance=None, created=False, **kwargs):
    """
    Delete the actual reference files when a reference info is deleted
    """
    try:
        os.unlink(f"{instance.file_location.path}.fxi")
    except FileNotFoundError as e:
        pass
    instance.file_location.delete(save=False)
    try:
        os.unlink(instance.minimap2_index_file_location)
    except FileNotFoundError as e:
        pass
