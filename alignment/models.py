from django.db import models

from reads.models import Barcode, FastqRead, FastqReadType, Run, GroupRun, Flowcell, JobMaster
from reference.models import ReferenceInfo, ReferenceLine


class PafRoughCov(models.Model):
    """
    Stores the PafCoverage in bins of ten, stores data for the Mapping Master/Detail chart
    """

    class Meta:
        indexes = [
            models.Index(fields=['bin_position_start'], name='position_index'),
        ]
    job_master = models.ForeignKey(
        JobMaster, on_delete=models.CASCADE, related_name="paf_rough_cov_list"
    )
    run = models.ForeignKey(
        Run, on_delete=models.CASCADE, related_name="prc_run", null=True, blank=True
    )
    flowcell = models.ForeignKey(
        Flowcell,
        on_delete=models.CASCADE,
        related_name="flowcell_prc_run",
        null=True,
        blank=True,
    )
    grouprun = models.ForeignKey(
        GroupRun,
        on_delete=models.CASCADE,
        related_name="pafroughcovs",
        null=True,
        blank=True,
    )
    read_type = models.ForeignKey(
        FastqReadType, related_name="prc_type", on_delete=models.SET_NULL, null=True
    )
    is_pass = models.BooleanField()  # pass = true, fail = false
    barcode = models.ForeignKey(
        Barcode, related_name="prc_barcode", null=True, on_delete=models.SET_NULL,
    )
    reference = models.ForeignKey(
        ReferenceInfo, related_name="paf_reference", on_delete=models.SET_NULL, null=True
    )
    chromosome = models.ForeignKey(
        ReferenceLine, related_name="paf_chromosome", on_delete=models.SET_NULL, null=True
    )
    bin_position_start = models.IntegerField(default=0)  # position
    bin_position_end = models.IntegerField(default=0)
    bin_coverage = models.IntegerField(default=0)
    # TODO I don't like this solution
    # TODO we could(should?) have this as a third table
    reference_pk = models.IntegerField(
        default=0
    )
    reference_name = models.CharField(
        max_length=256,
        null=True
    )
    chromosome_length = models.IntegerField(
        default=0
    )
    chromosome_pk = models.IntegerField(
        default=0
    )
    chromosome_name = models.CharField(
        max_length=256,
        null=True
    )

    def __str__(self):
        return f"{self.flowcell} {self.bin_position_start} - {self.bin_position_end}, Coverage: {self.bin_coverage}"


class PafSummaryCov(models.Model):
    """
    Store the aggregated results of the PafStore
    """

    job_master = models.ForeignKey(
        JobMaster, on_delete=models.CASCADE, related_name="paf_summary_cov_list"
    )
    barcode = models.ForeignKey(
        Barcode,
        on_delete=models.CASCADE,
        related_name="barcode_paf_summaries",
        null=True
    )
    barcode_name = models.CharField(max_length=32)
    chromosome = models.ForeignKey(
        ReferenceLine,
        related_name="paf_summary_chromosome",
        on_delete=models.SET_NULL,
        null=True,
    )
    # TODO to be deleted, use field from chromosome
    reference_line_length = models.IntegerField()
    read_count = models.BigIntegerField(default=0)
    total_yield = models.BigIntegerField(default=0)
    read_type = models.ForeignKey(
        FastqReadType,
        related_name="paf_summary_read_type",
        on_delete=models.SET_NULL,
        null=True,
    )
    coverage = models.FloatField(
        default=0
    )
    average_read_length = models.IntegerField(
        default=0
    )
    chromosome_pk = models.IntegerField(
        default=0
    )
    chromosome_name = models.CharField(
        max_length=256,
        null=True
    )
    reference_pk = models.IntegerField(
        default=0
    )
    reference_name = models.CharField(
        max_length=256,
        null=True
    )

    def __str__(self):
        if self.barcode:
            return f"Id: {self.id}, FlowcellId: {self.job_master.flowcell.id}, Barcode: {self.barcode.name}," \
                   f" Chromosome: {self.chromosome.line_name}"
        else:
            return f"Id: {self.id}, FlowcellId: {self.job_master.flowcell.id},"\
                   f" Chromosome: {self.chromosome.line_name}"
