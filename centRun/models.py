from django.db import models

# class MetaGenomicsMeta(models.Model):
#     """Meta data about the metagenomics run"""
#     task_name = models.CharField(max_length=200)
#     task_description = models.CharField(max_length=1000)
#     sample_description = models.CharField(max_length=1000)
#     timestamp = models.DateTimeField(default=timezone.now, null=True)
#     run_date = models.CharField(max_length=30, null=True)
#     run_time = models.CharField(max_length=30, null=True)
#     user = models.CharField(max_length=30, null=True)
#     # meta_id is the uuid used to link the centrifuge output sets with a particular job
#     meta_id = models.CharField(max_length=40, primary_key=True)
#     running = models.BooleanField()


class CentOutput(models.Model):
    """ Cent output to be put into the database. Fields include taxid, num reads, taskname, name, sumunique"""
    name = models.CharField(max_length=250, null=True)
    taxID = models.IntegerField(null=True)
    num_reads = models.IntegerField(null=True)
    sum_unique = models.IntegerField(null=True)
    taxa = models.CharField(null=True, max_length=100)
    # lineages = PickledObjectField(null=True)
    # uuid that links it to it's corresponding metadata
    task_meta = models.CharField(max_length=40)
    flowcell_id = models.IntegerField(null=True)
    # use djongo manager rather than django manager, in order to use mongo query syntax in run_centrifuge.py


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


class ReferenceGenomes(models.Model):
    """

    """
    species = models.CharField(max_length=65)
    refseq_link = models.CharField(max_length=100)


# tax_rank_filter = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
class LineageKey(models.Model):
    tax_id = models.IntegerField()


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
