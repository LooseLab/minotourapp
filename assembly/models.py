import datetime

from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.conf import settings

from reads.models import MinIONRun
from reads.models import FastqRead
from reads.models import FastqReadType
from reads.models import Barcode

# Create your models here.

class GfaStore(models.Model):
    run = models.ForeignKey(MinIONRun, on_delete=models.CASCADE, related_name='gfastorerun')
    readtype = models.ForeignKey(FastqReadType, related_name='readtypegfastore', blank=True,null=True)
    barcode = models.ForeignKey(Barcode,related_name='barcodegfastore', blank=True,null=True)
    nreads = models.IntegerField() # Number of input reads

    gfaformat = models.TextField()  #   string  The whole GFA file

    def __str__(self):
        return"{} {}".format(self.run,self.nreads)

class GfaSummary(models.Model):
    run = models.ForeignKey(MinIONRun, on_delete=models.CASCADE, related_name='gfasummaryrun')
    readtype = models.ForeignKey(FastqReadType, related_name='readtypegfasummary', blank=True,null=True)
    barcode = models.ForeignKey(Barcode,related_name='barcodegfasummary', blank=True,null=True)
    nreads = models.IntegerField() # Number of input reads
    ncontigs = models.IntegerField()#	int	Number of contigs
    maxlen = models.IntegerField()  #	int	Maximum contig length
    minlen  = models.IntegerField() #	int	Mininum contig length
    totlen = models.IntegerField()  #	int	Total contig length
    n50len = models.IntegerField() #    int Contig N50
    meanlen = models.IntegerField() #   int Mean contig length
    allcontigs = models.TextField() #   string array of contig lengths
    timecreated = models.DateTimeField(auto_now_add=True)

    def barcode_name(self):
        return self.barcode.name

    def type_name(self):
        return self.readtype.name

    def __str__(self):
        return"{} {} {}".format(self.run, self.barcode.name, self.nreads)
