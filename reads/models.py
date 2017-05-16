from django.db import models


class MinionRun(models.Model):
    run_name = models.CharField(max_length=64)
    run_id = models.CharField(max_length=64)
    barcode = models.CharField(max_length=32)

    def __str__(self):
        return self.run_name


class FastqReadType(models.Model):
    name = models.CharField(max_length=16)

    def __str__(self):
        return self.name


class FastqRead(models.Model):
    run_id = models.ForeignKey(MinionRun)
    read_id = models.CharField(max_length=64)
    read = models.IntegerField()
    channel = models.IntegerField()
    barcode = models.CharField(max_length=32)
    sequence = models.TextField()
    quality = models.TextField()
    is_pass = models.BooleanField() # pass = true, fail = false
    type = models.ForeignKey(FastqReadType)
    start_time = models.DateTimeField()

    def __str__(self):
        return self.read_id


class RunStatistic(models.Model):
    run_id = models.ForeignKey(MinionRun)
    sample_time = models.DateTimeField()
    total_length = models.IntegerField()
    max_length = models.IntegerField()
    min_length = models.IntegerField()
    average_length = models.DecimalField(max_digits=11, decimal_places=2)
    number_of_reads = models.IntegerField()
    number_of_channels = models.IntegerField()
    type = models.ForeignKey(FastqReadType)

    def __str__(self):
        return "{} {} {}".format(self.run_id,self.sample_time,self.type)

