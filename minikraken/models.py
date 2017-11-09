from django.db import models

# Create your models here.
from django.dispatch import receiver
from rest_framework.authtoken.models import Token
from reads.models import MinIONRun
from reads.models import FastqRead
from reads.models import FastqReadType
from reads.models import Barcode



class MiniKraken(models.Model):
    run = models.ForeignKey(MinIONRun, related_name='minikrakrun')
    read = models.ForeignKey(FastqRead, related_name='minikrakread')
    #reference = models.ForeignKey(ReferenceInfo, related_name='minikrakreference')
    #read_type = models.ForeignKey(FastqReadType, related_name='minikraktype',blank=True,null=True)
    #barcode = models.ForeignKey(Barcode, related_name='minikrakbarcode',blank=True,null=True)
    krakenstatus = models.CharField(max_length=1)
    krakentaxid = models.IntegerField()
    krakenseqlen = models.BigIntegerField()
    krakenlca = models.TextField()

    def __str__(self):
        return "{}".format(self.krakenlca)


class ParsedKraken(models.Model):
    run = models.ForeignKey(MinIONRun, related_name='parsedkrakrun')
    type = models.ForeignKey(FastqReadType, related_name='typekrakrun', blank=True,null=True)
    barcode = models.ForeignKey(Barcode,related_name='barcodekrakrun', blank=True,null=True)
    percentage = models.FloatField(null=True,blank=True)
    rootreads = models.IntegerField(null=True,blank=True)
    directreads = models.IntegerField(null=True,blank=True)
    rank = models.CharField(max_length=1,null=True,blank=True)
    NCBItaxid = models.IntegerField(null=True,blank=True)
    parent = models.CharField(max_length=128,null=True,blank=True)
    sci_name = models.CharField(max_length=128,null=True,blank=True)
    indentation = models.IntegerField(null=True,blank=True)
    orderin = models.IntegerField(null=True,blank=True)

    def barcode_name(self):
        return self.barcode.name

    def type_name(self):
        return self.type.name

    def __str__(self):
        return "{},{}".format(self.run,self.NCBItaxid)

