from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.conf import settings

from reads.models import MinIONRun
from reads.models import FastqRead
from reads.models import FastqReadType
from reads.models import Barcode

from reference.models import ReferenceInfo
from reference.models import ReferenceLine



# Create your models here.
##### Alignment Management


##### This class replicates original minoTour design but perhaps it is better just to store the sam file as a single row? Lets try that
"""
class SamStore(models.Model):
    run = models.ForeignKey(MinIONRun, related_name='runalignments')
    read = models.ForeignKey(FastqRead, related_name='readalignment')
    flag = models.IntegerField()
    rname = models.CharField(max_length=512)
    pos = models.BigIntegerField()
    mapq = models.IntegerField()
    cigar = models.TextField()
    rnext = models.CharField(max_length=512)
    pnext = models.IntegerField()
    tlen = models.IntegerField()
    seq = models.TextField()
    qual = models.TextField()
    nm = models.CharField(max_length=512)
    md = models.TextField()
    ass = models.CharField(max_length=512)
    xs = models.CharField(max_length=512)

    def __str__(self):
        return "{} {}".format(self.run_id,self.read_id)

"""

class PafStore(models.Model):
    run = models.ForeignKey(MinIONRun, related_name='pafalignemnts')
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

class SamStore(models.Model):
    run = models.ForeignKey(MinIONRun, related_name='runalignments')
    read = models.ForeignKey(FastqRead, related_name='readalignment')
    samline = models.TextField()

    def __str__(self):
        return "{} {}".format(self.run_id,self.read_id)

class SamRoughCov(models.Model):
    run = models.ForeignKey(MinIONRun, related_name='rc_runid')
    read_type = models.ForeignKey(FastqReadType, related_name='rc_type')
    #reference = models.TextField() #should switch to a reference database
    reference = models.ForeignKey(ReferenceInfo, related_name='samreference')
    #chromosome = models.TextField() #could be repeated of the reference name
    chromosome = models.ForeignKey(ReferenceLine, related_name='samchromosome')
    p = models.IntegerField() #position
    i = models.IntegerField() #incdel

    def __str__(self):
        return "{} {}".format(self.run,self.read_type)

class PafRoughCov(models.Model):
    run = models.ForeignKey(MinIONRun, related_name='prc_run')
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
    run = models.ForeignKey(MinIONRun, related_name='paf_summary')
    read_type = models.ForeignKey(FastqReadType, related_name='paf_summary_type')
    barcode = models.ForeignKey(Barcode, related_name='paf_summary_barcode', null=True)
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
            return (self.cumu_length/self.chrom_len())
        except AttributeError:
            return "undefined"

    def avg_read_len(self):
        try:
            return (self.cumu_length/self.read_count)
        except AttributeError:
            return "undefined"


    def __str__(self):
        return "{} {} {}".format(self.run,self.reference,self.chromosome)

@receiver(post_save, sender=PafStore)
def updatePafRoughCov(instance, sender, **kwargs):
    pafline = instance
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


@receiver(post_save, sender=SamStore)
def updateSamRoughCov(instance, sender, **kwargs):
    samline = instance
    print (samline)