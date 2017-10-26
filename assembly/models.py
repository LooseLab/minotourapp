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
    run = models.ForeignKey(MinIONRun, related_name='gfaassembly')
    nreads = models.IntegerField() # Number of input reads
    #GFA File Format summary:
    ncontigs = models.IntegerField()#	int	Number of contigs
    maxlen = models.IntegerField()  #	int	Maximum contig length
    minlen  = models.IntegerField() #	int	Mininum contig length
    totlen = models.IntegerField()  #	int	Total contig length
    n50len = models.IntegerField() #    int Contig N50
    meanlen = models.IntegerField() #   int Mean contig length
    gfaformat = models.TextField()  #   string  The whole GFA file

    def __str__(self):
        return"{} {}".format(self.run,self.nreads)
