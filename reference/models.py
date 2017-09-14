from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver

class ReferenceInfo(models.Model):
    reference_name = models.CharField(max_length=512)
    filename = models.CharField(max_length=512)
    bwa_index_file_location = models.CharField(max_length=512, null=True)
    minimap2_index_file_location = models.CharField(max_length=512, null=True)
    totalrefleN=models.IntegerField()

    def __str__(self):
        return "{} {}".format(self.reference_name, self.filename)

class ReferenceLine(models.Model):
    reference = models.ForeignKey(ReferenceInfo, related_name='referencelines')
    line_name = models.CharField(max_length=512)
    chromosome_length = models.IntegerField()

    def __str__(self):
        return "{} {}".format(self.reference,self.line_name)
