from django.db import models
from django.utils import timezone
from jobs.models import JobMaster


class MetaGenomicsMeta(models.Model):
    """"
    :purpose: Store information about the Metagenomics classification analysis
    :author: Rory

    Fields:

    :timestamp: (datetime) - Start of the analysis task
    :runtime: (str) - Time taken to run
    :flowcell_id: (int) - The Flowcell ID
    :task: (JobMaster Object) FK - the task record in the JobMaster that started this analysis
    :running: (bool) - If the analysis is currently running
    :number_of_reads: (int) - Number of reads in Flowcell
    :reads_classified: (int) - Number of reads that have been analysed
    :finish_time: (str) - The time the analysis finished

    """

    timestamp = models.DateTimeField(default=timezone.now, null=True)
    run_time = models.CharField(max_length=30, null=True)
    flowcell_id = models.IntegerField()
    task = models.ForeignKey(
        JobMaster,
        related_name="metadata_metagenomics",
    )
    running = models.BooleanField()
    number_of_reads = models.IntegerField()
    reads_classified = models.IntegerField()
    finish_time = models.CharField(max_length=40, null=True)


class CentOutput(models.Model):
    """
    :purpose: Store the summarised output of the metagenomics task, used to draw donut chart and All reads table
    :author: Rory
    Fields:

    :name: (str) - The name of the Species
    :tax_id: (int) - The taxonomic ID of this species
    :num_matches: (int) - The number of centrifuge matches to this species in this analysis
    :sum_unique: (int) - The number of centrifuge matches that match uniquely to this species
    :flowcell_id: (int) - The Id of the flowcell the reads came from #TODO make this a foreign key?
    :task: (JobMaster Object) FK - the task record in the JobMaster that started this analysis
    """
    name = models.CharField(max_length=250, null=True)
    tax_id = models.IntegerField(null=True)
    num_matches = models.IntegerField(null=True)
    sum_unique = models.IntegerField(null=True)
    flowcell_id = models.IntegerField(null=True)
    task = models.ForeignKey(
        JobMaster,
        related_name="centrifuge_summaries"
    )


class CartographyMapped(models.Model):
    """
        The traffic light species that have been identified
        # TODO currently ununsed
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
        # TODO currently unused
    """
    species = models.CharField(max_length=100)
    tax_id = models.IntegerField(null=True)
    cart_mapped_id = models.CharField(max_length=50)


class LineageValues(models.Model):
    """
        :purpose: Store the complete Lineages that are produced by NCBI taxa,
        manipulated into this more helpful format
        Used in centrifuge class, and all the views
        :author: Rory
        Fields:

        :tax_id: (int) - The taxonomic ID of this species
        :superkingdom: (str) - superkingdom member name if present
        :phylum: (str) - phylum member name if present
        :classy: (str) - class member name if present
        :order: (str) - order member name if present
        :family: (str) - fmaily member name if present
        :genus: (str) - genus member name if present
        :species: (str) - species member name if present
        :subspecies: (str) - subspecies member name if present
        :strain: (str) - strain member name if present
        :leaf: (str) - leaf member name if present
        :substrainspecies: (str) - subspecies strain member name if present

        """
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


class SankeyLinks(models.Model):
    """
        :purpose: Store the centrifuge classifier output in the format for the sankey diagram
        Used in centrifuge class, and all the sankey diagram view
        :author: Rory
        Fields:

        :source: (str) - The source node, that the link starts at, for example bacteria
        :target: (str) - The target node, that the link ends at, for example proteobacteria
        :value: (int) - The num of macthes illustrated by that link
        :tax_id: (int) - The tax_id of the target node
        :flowcell_id: (int) - The flowcell id that the centrifuger class was run on
        :task: (fk JobMAster) - The corresponding jobMaster object, used to seperate the datasets
    """
    source = models.CharField(null=True, max_length=100)
    target = models.CharField(null=True, max_length=100)
    value = models.IntegerField()
    tax_id = models.IntegerField()
    flowcell_id = models.IntegerField()
    task = models.ForeignKey(
        JobMaster,
        related_name="sankey_links"
    )
