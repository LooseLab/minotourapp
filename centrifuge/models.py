from django.db import models
from django.utils import timezone
from jobs.models import JobMaster


class MetaGenomicsMeta(models.Model):
    """Meta data about the metagenomics run"""

    timestamp = models.DateTimeField(default=timezone.now, null=True)
    run_time = models.CharField(max_length=30, null=True)
    flowcell_id = models.IntegerField()
    running = models.BooleanField()
    number_of_reads = models.IntegerField()
    reads_classified = models.IntegerField()
    finish_time = models.DateTimeField(null=True)


class CentOutput(models.Model):
    """ Cent output to be put into the database. Fields include taxid, num reads, taskname, name, sumunique"""
    name = models.CharField(max_length=250, null=True)
    tax_id = models.IntegerField(null=True)
    num_matches = models.IntegerField(null=True)
    sum_unique = models.IntegerField(null=True)
    taxa = models.CharField(null=True, max_length=100)
    flowcell_id = models.IntegerField(null=True)
    job_master = models.ForeignKey(
        JobMaster,
        related_name="centrifuge_summaries"
    )


class CartographyMapped(models.Model):
    """
        The traffic light species that have been identified
    """
    species = models.CharField(max_length=50)
    tax_id = models.IntegerField(null=True)
    total_reads = models.IntegerField(null=True)
    alert_level = models.IntegerField(null=True)
    num_reads = models.IntegerField(null=True)
    red_reads = models.IntegerField(null=True)
    task_meta = models.CharField(max_length=50)


class CartographyGuide(models.Model):
    """
        The species defaults for the detection
    """
    species = models.CharField(max_length=100)
    tax_id = models.IntegerField(null=True)
    cart_mapped_id = models.CharField(max_length=50)


class LineageValues(models.Model):
    tax_id = models.IntegerField()
    superkingdom = models.CharField(null=True, max_length=100)
    phylum = models.CharField(null=True, max_length=100)
    classy = models.CharField(null=True, max_length=100)
    order = models.CharField(null=True, max_length=100)
    family = models.CharField(null=True, max_length=100)
    genus = models.CharField(null=True, max_length=100)
    species = models.CharField(null=True, max_length=100)
    subspecies = models.CharField(null=True, max_length=100)
    strain = models.CharField(null=True, max_length=100)
    leaf = models.CharField(null=True, max_length=100)
    substrainspecies = models.CharField(null=True, max_length=100)
    task_meta = models.CharField(null=True, max_length=50)


class SankeyLinks(models.Model):
    source = models.CharField(null=True, max_length=100)
    target = models.CharField(null=True, max_length=100)
    value = models.IntegerField()
    tax_id = models.IntegerField()
    flowcell_id = models.IntegerField()
    rank = models.IntegerField()
