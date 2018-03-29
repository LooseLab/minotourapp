from django.db import models
# Create your models here.
from django.dispatch import receiver
from rest_framework.authtoken.models import Token

from devices.models import Flowcell
from reads.models import BarcodeGroup, FastqRead, FastqReadType, Run


class MiniKraken(models.Model):
    run = models.ForeignKey(Run, on_delete=models.CASCADE, related_name='minikrakrun', blank=True, null=True)
    flowcell = models.ForeignKey(Flowcell, on_delete=models.CASCADE, related_name='minikrakflowcell', blank=True, null=True)
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
    run = models.ForeignKey(Run, on_delete=models.CASCADE, related_name='parsedkrakrun', blank=True, null=True)
    flowcell = models.ForeignKey(Flowcell, on_delete=models.CASCADE, related_name='parsedkrakflowcell', blank=True, null=True)
    type = models.ForeignKey(FastqReadType, related_name='typekrakrun', blank=True,null=True)
    barcode = models.ForeignKey(BarcodeGroup,related_name='barcodekrakrun', blank=True,null=True)
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
        if self.run is not None:
            return "{} {}".format(self.run, self.NCBItaxid)
        else:
            #return "{} {}".format(self.flowcell, self.NCBItaxid)
            return "{} {}".format(self.flowcell_id, self.NCBItaxid)
