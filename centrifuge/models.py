"""models.py

"""
from django.db import models
from django.utils import timezone
from jobs.models import JobMaster
from reads.models import Barcode
from reads.models import Flowcell


class SankeyLink(models.Model):
    """
        :purpose: Store the centrifuge classifier output in the format for the sankey diagram
        Used in centrifuge class, and all the sankey diagram view
        :author: Rory
        Fields:

        :source: (str) - The source node, that the link starts at, for example bacteria
        :target: (str) - The target node, that the link ends at, for example proteobacteria
        :value: (int) - The num of matches illustrated by that link
        :tax_id: (int) - The tax_id of the target node
        :flowcell_id: (int) - The flowcell id that the centrifuger class was run on
        :task: (fk JobMaster) - The corresponding jobMaster object, used to separate the datasets
        :target_tax_level: - The taxa level of the target node, i.e Phylum for Bacteria (kingdom) \
        to ProteoBacteria (Phylum)
    """
    source = models.CharField(null=True, max_length=100)

    target = models.CharField(null=True, max_length=100)

    value = models.IntegerField(null=True, default=0)

    barcode_name = models.CharField(null=True, max_length=100)

    tax_id = models.IntegerField()

    task = models.ForeignKey(
        JobMaster,
        related_name="sankey_link"
    )

    target_tax_level = models.CharField(max_length=100)

    path = models.IntegerField(default=0)

    def __str__(self):
        return "{} {} {} {}".format(self.barcode_name, self.tax_id, self.target_tax_level, self.task)


class Metadata(models.Model):
    """"
        :purpose: Store information about the Metagenomics classification analysis, used in centrifuger.py
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

    task = models.ForeignKey(
        JobMaster,
        related_name="metadata",
    )

    start_time = models.DateTimeField(default=timezone.now, null=True)

    finish_time = models.CharField(max_length=40, null=True)


class CentrifugeOutput(models.Model):
    """
        :purpose: Store the unsummarised output of the metagenomics task for all reads
        ,to be parsed and put into the next table
        :author: Rory
        Fields:

        :name: (str) - The name of the Species
        :tax_id: (int) - The taxonomic ID of this species
        :num_matches: (int) - The number of centrifuge matches to this species in this analysis
        :sum_unique: (int) - The number of centrifuge matches that match uniquely to this species
        :flowcell_id: (int) - The Id of the flowcell the reads came from
        :task: (JobMaster Object) FK - the task record in the JobMaster that started this analysis
    """

    name = models.CharField(max_length=250, null=True)

    tax_id = models.IntegerField(null=True)

    task = models.ForeignKey(
        JobMaster,
        related_name="centrifuge_output",
        on_delete=models.CASCADE,
    )

    num_matches = models.IntegerField(default=0)

    sum_unique = models.IntegerField(default=0)

    barcode_name = models.CharField(max_length=50)

    proportion_of_classified = models.FloatField(default=0.0)

    superkingdom = models.CharField(null=True, max_length=70)

    phylum = models.CharField(null=True, max_length=70)

    classy = models.CharField(null=True, max_length=70)

    order = models.CharField(null=True, max_length=70)

    family = models.CharField(null=True, max_length=70)

    genus = models.CharField(null=True, max_length=70)

    species = models.CharField(null=True, max_length=70)

    latest = models.IntegerField(null=True)

    def __str__(self):
        return "{} {} {} {}".format(self.barcode_name, self.tax_id, self.name, self.task)


class CentrifugeOutputParsed(models.Model):
    """
        :purpose: Store the summarised output of the metagenomics task for all reads
        , used to draw donut chart and All reads table
        :author: Rory
        Fields:

        :name: (str) - The name of the Species
        :tax_id: (int) - The taxonomic ID of this species
        :num_matches: (int) - The number of centrifuge matches to this species in this analysis
        :sum_unique: (int) - The number of centrifuge matches that match uniquely to this species
        :flowcell_id: (int) - The Id of the flowcell the reads came from
        :task: (JobMaster Object) FK - the task record in the JobMaster that started this analysis
    """

    name = models.CharField(max_length=250, null=True)

    tax_id = models.IntegerField(null=True)

    task = models.ForeignKey(
        JobMaster,
        related_name="centrifuge_output_parsed",
        on_delete=models.CASCADE,
    )

    num_matches = models.IntegerField(default=0)

    sum_unique = models.IntegerField(default=0)

    barcode_name = models.CharField(max_length=50)

    proportion_of_classified = models.FloatField(default=0.0)

    superkingdom = models.CharField(null=True, max_length=70)

    phylum = models.CharField(null=True, max_length=70)

    classy = models.CharField(null=True, max_length=70)

    order = models.CharField(null=True, max_length=70)

    family = models.CharField(null=True, max_length=70)

    genus = models.CharField(null=True, max_length=70)

    species = models.CharField(null=True, max_length=70)

    latest = models.IntegerField(default=0)

    def __str__(self):
        return "{} {} {} {}".format(self.barcode_name, self.tax_id, self.name, self.task)


class MappingResult(models.Model):
    """
        The traffic light species that have been identified
    """
    task = models.ForeignKey(
        JobMaster,
        related_name="mapping_result",
        on_delete=models.CASCADE,
        null=True
    )

    species = models.CharField(max_length=50)

    tax_id = models.IntegerField(null=True)

    barcode_name = models.CharField(max_length=50, default="All reads")

    num_mapped = models.IntegerField(default=0)

    mapped_proportion_of_classified = models.FloatField(default=0)

    red_reads = models.IntegerField(null=True, default=0)

    red_reads_proportion_of_classified = models.FloatField(default=0)

    num_matches = models.IntegerField(default=0)

    proportion_of_classified = models.FloatField(default=0.0)

    sum_unique = models.IntegerField(default=0)

    def __str__(self):
        return "{} {} {}".format(self.species, self.tax_id, self.task)


class TargetMappedReadId(models.Model):
    """
        The read ids for reads that have identified dangerously
    """
    read_id = models.CharField(max_length=100, unique=True)

    mapping_result = models.ForeignKey(
        MappingResult,
        related_name="target_read_id",
        on_delete=models.CASCADE
    )

    def __str__(self):
        return self.read_id


class MappingTarget(models.Model):
    """
        The species defaults for the detection
    """
    species = models.CharField(max_length=100)

    tax_id = models.IntegerField(null=True)

    target_set = models.CharField(default="Default", max_length=100)

    start = models.IntegerField(default=0)

    end = models.IntegerField(default=1)

    name = models.CharField(default="danger_zone", max_length=50)

    gff_line_type = models.CharField(default="gene", max_length=100)

    def __str__(self):
        return "{} {} {}-{} {}".format(self.species, self.target_set, self.start, self.end, self.gff_line_type)


class LineageValue(models.Model):
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


class DonutData(models.Model):
    """
        The data for the donut chart and the donut ranking table
    """
    task = models.ForeignKey(
        JobMaster,
        on_delete=models.CASCADE,
        related_name="donut_data"
    )
    num_matches = models.IntegerField(default=0)
    sum_unique = models.IntegerField(default=0)
    name = models.CharField(max_length=60)
    barcode_name = models.CharField(max_length=40, default="All reads")
    tax_rank = models.CharField(max_length=80, default="Species")
    latest = models.IntegerField(null=True, default=0)

    def __str__(self):
        return "{} {} {}".format(self.task, self.name, self.barcode_name)
