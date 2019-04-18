from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver

from jobs.models import JobMaster
from reads.models import Barcode, FastqRead, FastqReadType, Run, GroupRun, Flowcell
from reference.models import ReferenceInfo, ReferenceLine


class PafStore(models.Model):

    job_master = models.ForeignKey(

        JobMaster,
        on_delete=models.CASCADE,
        related_name='pafstore_list'
    )

    flowcell = models.ForeignKey(  # TODO dele - read from job_master.flowcell

        Flowcell,
        on_delete=models.CASCADE,
        related_name='flowcellpafalignemnts',
        null=True,
        blank=True
    )

    read = models.ForeignKey(

        FastqRead,
        related_name='pafreadalignment',
        on_delete=models.DO_NOTHING,
    )

    reference = models.ForeignKey(  # TODO delete - read from job_master.reference

        ReferenceInfo,
        related_name='pafstorereference',
        on_delete=models.DO_NOTHING,
    )

    # PAF File Format:

    qsn = models.CharField(

        max_length=256
    )  # 1	string	Query sequence name

    qsl = models.IntegerField(

    )  # 2	int	Query sequence length

    qs = models.IntegerField(

    )  # 3	int	Query start (0-based)

    qe = models.IntegerField(

    )  # 4	int	Query end (0-based)

    rs = models.CharField(

        max_length=1
    )  # 5	char	Relative strand: "+" or "-"

    tsn = models.ForeignKey(

        ReferenceLine,
        on_delete=models.CASCADE,
        related_name='pafstore_list'
    )  # 6	string	Target sequence name

    tsl = models.IntegerField(

    )  # 7	int	Target sequence length

    ts = models.IntegerField(

    )  # 8	int	Target start on original strand (0-based)

    te = models.IntegerField(

    )  # 9	int	Target end on original strand (0-based)

    nrm = models.IntegerField(

    )  # 10	int	Number of residue matches

    abl = models.IntegerField(

    )  # 11	int	Alignment block length

    mq = models.IntegerField(

    )  # 12	int	Mapping quality (0-255; 255 for missing)

    def __str__(self):
        return "{}".format(self.qsn)


class PafRoughCov(models.Model):

    job_master = models.ForeignKey(

        JobMaster,
        on_delete=models.CASCADE,
        related_name='pafroughcov_list'
    )

    run = models.ForeignKey(
        Run,
        on_delete=models.CASCADE,
        related_name='prc_run',
        null=True,
        blank=True
    )

    flowcell = models.ForeignKey(
        Flowcell,
        on_delete=models.CASCADE,
        related_name='flowcell_prc_run',
        null=True,
        blank=True
    )

    grouprun = models.ForeignKey(
        GroupRun,
        on_delete=models.CASCADE,
        related_name='pafroughcovs',
        null=True,
        blank=True
    )

    read_type = models.ForeignKey(
        FastqReadType,
        related_name='prc_type',
        on_delete=models.DO_NOTHING,
    )

    is_pass = models.BooleanField(

    )  # pass = true, fail = false

    barcode = models.ForeignKey(
        Barcode,
        related_name='prc_barcode',
        null=True,
        on_delete=models.DO_NOTHING,
    )

    barcode_name = models.CharField(

        max_length=32
    )

    reference = models.ForeignKey(
        ReferenceInfo,
        related_name='pafreference',
        on_delete=models.DO_NOTHING,
    )

    chromosome = models.ForeignKey(
        ReferenceLine,
        related_name='pafchromosome',
        on_delete=models.DO_NOTHING,
    )

    p = models.IntegerField(

    )  # position

    i = models.IntegerField(
        default=0
    )  # incdel

    def to_chart_data(self):

        return [ self.p, self.i]

    def __str__(self):
        return "{} {}".format(self.run, self.p)


class PafSummaryCov(models.Model):

    job_master = models.ForeignKey(

        JobMaster,
        on_delete=models.CASCADE,
        related_name='paf_summary_cov_list'
    )

    barcode_name = models.CharField(

        max_length=32
    )

    reference_line_name = models.CharField(

        max_length=512
    )

    reference_line_length = models.BigIntegerField(

    )

    read_count = models.BigIntegerField(

        default=0
    )

    total_length = models.BigIntegerField(

        default=0
    )


@receiver(post_save, sender=PafStore)
def updatePafRoughCov(instance, sender, **kwargs):
    pafline = instance
    if pafline.run is not None:
        pafstart, created1 = PafRoughCov.objects.update_or_create(
            run=pafline.run,
            read_type=pafline.read.type,
            barcode=pafline.read.barcode,
            reference=pafline.reference,
            chromosome=pafline.tsn,
            p=pafline.ts  # position
        )
        pafstart.i += 1
        pafstart.save()
        pafend, created2 = PafRoughCov.objects.update_or_create(
            run=pafline.run,
            read_type=pafline.read.type,
            barcode=pafline.read.barcode,
            reference=pafline.reference,
            chromosome=pafline.tsn,
            p=(pafline.te) + 1  # position
        )
        pafend.i -= 1
        pafend.save()
    else:
        pafstart, created1 = PafRoughCov.objects.update_or_create(
            flowcell=pafline.flowcell,
            read_type=pafline.read.type,
            barcode=pafline.read.barcode,
            reference=pafline.reference,
            chromosome=pafline.tsn,
            p=pafline.ts  # position
        )
        pafstart.i += 1
        pafstart.save()
        pafend, created2 = PafRoughCov.objects.update_or_create(
            flowcell=pafline.flowcell,
            read_type=pafline.read.type,
            barcode=pafline.read.barcode,
            reference=pafline.reference,
            chromosome=pafline.tsn,
            p=(pafline.te) + 1  # position
        )
        pafend.i -= 1
        pafend.save()
