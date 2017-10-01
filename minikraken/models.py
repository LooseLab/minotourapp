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
    #read_type = models.ForeignKey(FastqReadType, related_name='minikraktype')
    krakenstatus = models.CharField(max_length=1)
    krakentaxid = models.IntegerField()
    krakenseqlen = models.BigIntegerField()
    krakenlca = models.TextField()

    def __str__(self):
        return "{}".format(self.krakenlca)


class TestKraken(models.Model):
    stuff1 = models.CharField(max_length=20)


    def __str__(self):
        return "{}".format(self.stuff1)

