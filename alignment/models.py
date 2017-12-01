from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver

from reads.models import Barcode
from reads.models import BarcodeGroup
from reads.models import FastqRead
from reads.models import FastqReadType
from reads.models import MinIONRun
from reads.models import FlowCell
from reference.models import ReferenceInfo
from reference.models import ReferenceLine


class PafStore(models.Model):
    run = models.ForeignKey(MinIONRun, on_delete=models.CASCADE, related_name='pafalignemnts', null=True, blank=True)
    flowcell = models.ForeignKey(FlowCell, on_delete=models.CASCADE, related_name='flowcellpafalignemnts', null=True,blank=True)
    read = models.ForeignKey(FastqRead, related_name='pafreadalignment')
    reference = models.ForeignKey(ReferenceInfo, related_name='pafstorereference')
    read_type = models.ForeignKey(FastqReadType, related_name='pafstoretype')
    #PAF File Format:
    qsn = models.CharField(max_length=256)#1	string	Query sequence name
    qsl = models.IntegerField()#2	int	Query sequence length
    qs  = models.IntegerField()#3	int	Query start (0-based)
    qe = models.IntegerField()#4	int	Query end (0-based)
    rs = models.CharField(max_length=1)#5	char	Relative strand: "+" or "-"
    #tsn = models.CharField(max_length=256)#6	string	Target sequence name
    tsn = models.ForeignKey(ReferenceLine, related_name='pafstorechromosome')#6	string	Target sequence name
    tsl = models.IntegerField()#7	int	Target sequence length
    ts = models.IntegerField()#8	int	Target start on original strand (0-based)
    te = models.IntegerField()#9	int	Target end on original strand (0-based)
    nrm = models.IntegerField()#10	int	Number of residue matches
    abl = models.IntegerField()#11	int	Alignment block length
    mq = models.IntegerField()#12	int	Mapping quality (0-255; 255 for missing)

    def __str__(self):
        return"{} {}".format(self.run,self.qsn)


class PafStore_transcriptome(models.Model):
    run = models.ForeignKey(MinIONRun, on_delete=models.CASCADE, related_name='pafalignemnts_transcriptome')
    read = models.ForeignKey(FastqRead, related_name='pafreadalignment_transcriptome')
    reference = models.ForeignKey(ReferenceInfo, related_name='pafstorereference_transcriptome')
    read_type = models.ForeignKey(FastqReadType, related_name='pafstoretype_transcriptome')
    #PAF File Format:
    qsn = models.CharField(max_length=256)#1	string	Query sequence name
    qsl = models.IntegerField()#2	int	Query sequence length
    qs  = models.IntegerField()#3	int	Query start (0-based)
    qe = models.IntegerField()#4	int	Query end (0-based)
    rs = models.CharField(max_length=1)#5	char	Relative strand: "+" or "-"
    #tsn = models.CharField(max_length=256)#6	string	Target sequence name
    tsn = models.ForeignKey(ReferenceLine, related_name='pafstorechromosome_transcriptome')#6	string	Target sequence name
    tsl = models.IntegerField()#7	int	Target sequence length
    ts = models.IntegerField()#8	int	Target start on original strand (0-based)
    te = models.IntegerField()#9	int	Target end on original strand (0-based)
    nrm = models.IntegerField()#10	int	Number of residue matches
    abl = models.IntegerField()#11	int	Alignment block length
    mq = models.IntegerField()#12	int	Mapping quality (0-255; 255 for missing)

    def __str__(self):
        return"{} {}".format(self.run,self.qsn)

class SamStore(models.Model):
    run = models.ForeignKey(
        MinIONRun,
        related_name='runalignments'
    )

    read = models.ForeignKey(
        FastqRead,
        related_name='readalignment'
    )

    samline = models.TextField(

    )

    def __str__(self):
        return "{} {}".format(self.run_id, self.read_id)


class SamRoughCov(models.Model):
    run = models.ForeignKey(
        MinIONRun,
        related_name='rc_runid'
    )

    read_type = models.ForeignKey(
        FastqReadType,
        related_name='rc_type'
    )

    reference = models.ForeignKey(
        ReferenceInfo,
        related_name='samreference'
    )

    chromosome = models.ForeignKey(
        ReferenceLine,
        related_name='samchromosome'
    )

    p = models.IntegerField(

    ) #position

    i = models.IntegerField(

    ) #incdel

    def __str__(self):
        return "{} {}".format(self.run,self.read_type)


class PafRoughCov(models.Model):
    run = models.ForeignKey(MinIONRun, on_delete=models.CASCADE, related_name='prc_run', null=True,
                                 blank=True)
    flowcell = models.ForeignKey(FlowCell, on_delete=models.CASCADE, related_name='flowcell_prc_run', null=True,
                                 blank=True)
    read_type = models.ForeignKey(FastqReadType, related_name='prc_type')
    barcode = models.ForeignKey(Barcode, related_name='prc_barcode', null=True)
    #reference = models.TextField() #should switch to a reference database
    reference = models.ForeignKey(ReferenceInfo, related_name='pafreference')
    #chromosome = models.TextField() #could be repeat of the reference name
    chromosome = models.ForeignKey(ReferenceLine, related_name='pafchromosome')
    p = models.IntegerField() #position
    i = models.IntegerField(default=0) #incdel

    def __str__(self):
        return "{} {}".format(self.run,self.p)


class PafSummaryCov(models.Model):
    run = models.ForeignKey(MinIONRun, on_delete=models.CASCADE, related_name='paf_summary', null=True,blank=True)
    flowcell = models.ForeignKey(FlowCell, on_delete=models.CASCADE, related_name='flowcell_paf_summary', null=True,blank=True)
    read_type = models.ForeignKey(FastqReadType, related_name='paf_summary_type')
    barcode = models.ForeignKey(Barcode, related_name='paf_summary_barcode', null=True)
    barcodegroup = models.ForeignKey(BarcodeGroup, related_name='paf_summary_barcodegroup', null=True)
    reference = models.ForeignKey(ReferenceInfo,related_name='paf_summary_reference')
    chromosome = models.ForeignKey(ReferenceLine, related_name='paf_summary_chromosome')
    read_count = models.BigIntegerField(default=0)
    cumu_length = models.BigIntegerField(default=0)

    def barcode_name(self):
        try:
            return self.barcode.name
        except AttributeError:
            return "undefined"

    def read_type_name(self):
        try:
            return self.read_type.name
        except AttributeError:
            return "undefined"

    def ref_name(self):
        try:
            return self.reference.reference_name
        except AttributeError:
            return "undefined"

    def ref_len(self):
        try:
            return self.reference.totalrefleN
        except AttributeError:
            return "undefined"

    def chrom_name(self):
        try:
            return self.chromosome.line_name
        except AttributeError:
            return "undefined"

    def chrom_len(self):
        try:
            return self.chromosome.chromosome_length
        except AttributeError:
            return "undefined"

    def chrom_cover(self):
        try:
            return (round(self.cumu_length/self.chrom_len(),3))
        except AttributeError:
            return "undefined"

    def avg_read_len(self):
        try:
            return (round(self.cumu_length/self.read_count,3))
        except AttributeError:
            return "undefined"


    def __str__(self):
        return "{} {} {}".format(self.run,self.reference,self.chromosome)


class PafRoughCov_transcriptome(models.Model):
    run = models.ForeignKey(MinIONRun, on_delete=models.CASCADE, related_name='prc_run_transcriptome')
    read_type = models.ForeignKey(FastqReadType, related_name='prc_type_transcriptome')
    barcode = models.ForeignKey(Barcode, related_name='prc_barcode_transcriptome', null=True)
    #reference = models.TextField() #should switch to a reference database
    reference = models.ForeignKey(ReferenceInfo, related_name='pafreference_transcriptome')
    #chromosome = models.TextField() #could be repeat of the reference name
    chromosome = models.ForeignKey(ReferenceLine, related_name='pafchromosome_transcriptome')
    p = models.IntegerField() #position
    i = models.IntegerField(default=0) #incdel

    def __str__(self):
        return "{} {}".format(self.run,self.p)

class PafSummaryCov_transcriptome(models.Model):
    run = models.ForeignKey(MinIONRun, on_delete=models.CASCADE, related_name='paf_summary_transcriptome')
    read_type = models.ForeignKey(FastqReadType, related_name='paf_summary_type_transcriptome')
    barcode = models.ForeignKey(Barcode, related_name='paf_summary_barcode_transcriptome', null=True)
    reference = models.ForeignKey(ReferenceInfo,related_name='paf_summary_reference_transcriptome')
    chromosome = models.ForeignKey(ReferenceLine, related_name='paf_summary_chromosome_transcriptome')
    read_count = models.BigIntegerField(default=0)
    cumu_length = models.BigIntegerField(default=0)

    def barcode_name(self):
        try:
            return self.barcode.name
        except AttributeError:
            return "undefined"

    def read_type_name(self):
        try:
            return self.read_type.name
        except AttributeError:
            return "undefined"

    def ref_name(self):
        try:
            return self.reference.reference_name
        except AttributeError:
            return "undefined"

    def ref_len(self):
        try:
            return self.reference.totalrefleN
        except AttributeError:
            return "undefined"

    def chrom_name(self):
        try:
            return self.chromosome.line_name
        except AttributeError:
            return "undefined"

    def chrom_len(self):
        try:
            return self.chromosome.chromosome_length
        except AttributeError:
            return "undefined"

    def chrom_cover(self):
        try:
            return (round(self.cumu_length/self.chrom_len(),3))
        except AttributeError:
            return "undefined"

    def avg_read_len(self):
        try:
            return (round(self.cumu_length/self.read_count,3))
        except AttributeError:
            return "undefined"


    def __str__(self):
        return "{} {} {}".format(self.run,self.reference,self.chromosome)


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
            p=pafline.ts #position
        )
        pafstart.i+=1
        pafstart.save()
        pafend, created2 = PafRoughCov.objects.update_or_create(
            run=pafline.run,
            read_type=pafline.read.type,
            barcode=pafline.read.barcode,
            reference=pafline.reference,
            chromosome=pafline.tsn,
            p=(pafline.te)+1 #position
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


@receiver(post_save, sender=PafStore_transcriptome)
def updatePafRoughCov_transcriptome(instance, sender, **kwargs):
    pafline = instance
    pafstart, created1 = PafRoughCov_transcriptome.objects.update_or_create(
        run=pafline.run,
        read_type=pafline.read.type,
        barcode=pafline.read.barcode,
        reference=pafline.reference,
        chromosome=pafline.tsn,
        p=pafline.ts #position
    )
    pafstart.i+=1
    pafstart.save()
    pafend, created2 = PafRoughCov_transcriptome.objects.update_or_create(
        run=pafline.run,
        read_type=pafline.read.type,
        barcode=pafline.read.barcode,
        reference=pafline.reference,
        chromosome=pafline.tsn,
        p=(pafline.te)+1 #position
    )
    pafend.i -= 1
    pafend.save()


@receiver(post_save, sender=SamStore)
def updateSamRoughCov(instance, sender, **kwargs):
    samline = instance
    print (samline)
