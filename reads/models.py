import datetime

from django.db import models
from django.conf import settings
from django.db.models.signals import post_save
from django.dispatch import receiver
from rest_framework.authtoken.models import Token


class MinION(models.Model):
    minION_name = models.CharField(max_length=64)
    owner = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='minIONs')

    def __str__(self):
        return self.minION_name

    def status(self):
        return self.events.order_by('datetime').last()


class MinIONRun(models.Model):
    run_name = models.CharField(max_length=64)
    run_id = models.CharField(max_length=64)
    is_barcoded = models.BooleanField(default=False)
    owner = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='runs')
    minION = models.ForeignKey(MinION, blank=True, null=True)

    def __str__(self):
        return self.run_name


class MinIONEventType(models.Model):
    name = models.CharField(max_length=64)

    def __str__(self):
        return self.name


class MinIONEvent(models.Model):
    computer_name = models.CharField(max_length=256)
    minION = models.ForeignKey(MinION, related_name='events')
    event = models.ForeignKey(MinIONEventType)
    datetime =models.DateTimeField()

    def __str__(self):
        return "{} {} {} {}".format(self.computer_name, self.minION, self.event, self.datetime)


class FastqReadType(models.Model):
    name = models.CharField(max_length=16)

    def __str__(self):
        return self.name


class FastqRead(models.Model):
    run_id = models.ForeignKey(MinIONRun,
                               on_delete=models.CASCADE,
                               related_name='reads')
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


class RunSummary(models.Model):
    run_id = models.ForeignKey(MinIONRun,on_delete=models.CASCADE)
    total_length = models.IntegerField(default=0)
    read_count = models.IntegerField(default=0)
    type = models.ForeignKey(FastqReadType)
    max_length = models.IntegerField(default=0)
    min_length = models.IntegerField(default=0)
    pass_length = models.IntegerField(default=0)
    pass_max_length = models.IntegerField(default=0)
    pass_min_length = models.IntegerField(default=0)
    pass_count = models.IntegerField(default=0)

    def __str__(self):
        return "{} {} {} {}".format(self.run_id, self.total_length, self.read_count, self.type)


class RunSummaryBarCode(models.Model):
    run_id = models.ForeignKey(MinIONRun,on_delete=models.CASCADE)
    total_length = models.IntegerField(default=0)
    read_count = models.IntegerField(default=0)
    type = models.ForeignKey(FastqReadType)
    barcode = models.CharField(max_length=32)
    pass_length = models.IntegerField(default=0)
    pass_max_length = models.IntegerField(default=0)
    pass_min_length = models.IntegerField(default=0)
    pass_count = models.IntegerField(default=0)
    max_length = models.IntegerField(default=0)
    min_length = models.IntegerField(default=0)

    def __str__(self):
        return "{} {} {} {} {}".format(self.run_id, self.total_length, self.read_count, self.type, self.barcode)


class RunStatistic(models.Model):
    run_id = models.ForeignKey(MinIONRun, on_delete=models.CASCADE)
    sample_time = models.DateTimeField()
    total_length = models.IntegerField(default=0)
    max_length = models.IntegerField(default=0)
    min_length = models.IntegerField(default=0)
    #average_length = models.DecimalField(max_digits=11, decimal_places=2)
    read_count = models.IntegerField(default=0)
    pass_length = models.IntegerField(default=0)
    pass_max_length = models.IntegerField(default=0)
    pass_min_length = models.IntegerField(default=0)
    pass_count = models.IntegerField(default=0)
    #number_of_channels = models.IntegerField() #Need to work out how to implement this!
    type = models.ForeignKey(FastqReadType)

    def __str__(self):
        return "{} {} {}".format(self.run_id, self.sample_time, self.type)


class RunStatisticBarcode(models.Model):
    run_id = models.ForeignKey(MinIONRun, on_delete=models.CASCADE)
    sample_time = models.DateTimeField()
    total_length = models.IntegerField(default=0)
    read_count = models.IntegerField(default=0)
    max_length = models.IntegerField(default=0)
    min_length = models.IntegerField(default=0)
    #average_length = models.DecimalField(default=0,max_digits=11, decimal_places=2)
    pass_length = models.IntegerField(default=0)
    pass_max_length = models.IntegerField(default=0)
    pass_min_length = models.IntegerField(default=0)
    pass_count = models.IntegerField(default=0)
    #number_of_channels = models.IntegerField() #Need to work out how to implement this!
    type = models.ForeignKey(FastqReadType)
    barcode = models.CharField(max_length=32)

    def __str__(self):
        return "{} {} {} {}".format(self.run_id, self.sample_time, self.type, self.barcode)


@receiver(post_save, sender=settings.AUTH_USER_MODEL)
def create_auth_token(sender, instance=None, created=False, **kwargs):
    if created:
        Token.objects.create(user=instance)



@receiver(post_save, sender=FastqRead)
def update_global_state(instance, sender, **kwargs):
    ipn_obj = instance
    #Here we work out the start time
    tm = ipn_obj.start_time

    tm = tm - datetime.timedelta(minutes=(tm.minute % 1) - 1,
                                 seconds=tm.second,
                                 microseconds=tm.microsecond)
    #print (type(tm))

    obj1, created = RunSummary.objects.update_or_create(
        run_id=ipn_obj.run_id, type=ipn_obj.type #, defaults={'total_length': 0, 'read_count': 0}
    )
    update_sum_stats(obj1, ipn_obj)

    obj2,created2 = RunSummaryBarCode.objects.update_or_create(
        run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=ipn_obj.barcode  # , defaults={'total_length': 0, 'read_count': 0}
    )
    update_sum_stats(obj2, ipn_obj)

    obj3, created3 = RunStatistic.objects.update_or_create(
        run_id=ipn_obj.run_id, type=ipn_obj.type, sample_time=tm
    )
    update_sum_stats(obj3, ipn_obj)

    obj4, created4 = RunStatisticBarcode.objects.update_or_create(
        run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=ipn_obj.barcode, sample_time=tm
        # , defaults={'total_length': 0, 'read_count': 0}
    )
    update_sum_stats(obj4, ipn_obj)


def update_sum_stats(obj,ipn_obj):
    if ipn_obj.is_pass:
        obj.pass_length += len(ipn_obj.sequence)
        if len(ipn_obj.sequence) > obj.pass_max_length:
            obj.pass_max_length = len(ipn_obj.sequence)
        if obj.pass_min_length == 0:
            obj.pass_min_length = len(ipn_obj.sequence)
        if len(ipn_obj.sequence) < obj.pass_min_length:
            obj.pass_min_length = len(ipn_obj.sequence)
        obj.pass_count += 1
    obj.total_length += len(ipn_obj.sequence)
    if len(ipn_obj.sequence) > obj.max_length:
        obj.max_length = len(ipn_obj.sequence)
    if obj.min_length == 0:
        obj.min_length = len(ipn_obj.sequence)
    if len(ipn_obj.sequence) < obj.min_length:
        obj.min_length = len(ipn_obj.sequence)
    obj.read_count += 1
    obj.save()
