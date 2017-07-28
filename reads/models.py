import datetime

from django.conf import settings
from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver
from rest_framework.authtoken.models import Token


class MinION(models.Model):
    minION_name = models.CharField(max_length=64)
    owner = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='minIONs')

    class Meta:
        verbose_name = 'MinION'
        verbose_name_plural = 'MinIONs'

    def __str__(self):
        return self.minION_name

    def status(self):
        try:
            status = self.events.order_by('datetime').last().event.name
        except AttributeError:
            status = "undefined"
        return status

    def last_run(self):
        try:
            return self.minionrun.last().id
            #return reverse('minIONrunstats_list', args=[self.minionrun.last().id])
        except AttributeError:
            last_run = "undefined"
        return last_run

    def computer(self):
        try:
            computer = self.events.order_by('datetime').last().computer_name
        except AttributeError:
            computer = "undefined"
        return computer

    def sample_name(self):
        try:
            return self.currentdetails.minKNOW_sample_name
        except AttributeError:
            return "undefined"

    def minKNOW_version(self):
        try:
            return self.currentrundetails.minKNOW_version
        except AttributeError:
            return "undefined"

    def flow_cell_id(self):
        try:
            return self.currentdetails.minKNOW_flow_cell_id
        except AttributeError:
            return "undefined"

    def run_status(self):
        try:
            return self.currentdetails.minKNOW_status
        except AttributeError:
            return "undefined"

    def run_name(self):
        try:
            return self.currentdetails.minKNOW_run_name
        except AttributeError:
            return "undefined"

    def total_drive_space(self):
        try:
            return self.currentdetails.minKNOW_total_drive_space
        except AttributeError:
            return "undefined"

    def space_till_shutdown(self):
        try:
            return self.currentdetails.minKNOW_disk_space_till_shutdown
        except AttributeError:
            return "undefined"

    def space_available(self):
        try:
            return self.currentdetails.minKNOW_disk_available
        except AttributeError:
            return "undefined"

    def warnings(self):
        try:
            return self.currentdetails.minKNOW_warnings
        except AttributeError:
            return "undefined"

    def currentscript(self):
        try:
            return self.currentdetails.minKNOW_current_script
        except AttributeError:
            return "undefined"

    def event_yield(self):
        try:
            return self.currentrunstats.order_by('sample_time').last().event_yield
        except AttributeError:
            return 0

    def voltage_value(self):
        try:
            return self.currentrunstats.order_by('sample_time').last().voltage_value
        except AttributeError:
            return 0


class MinIONControl(models.Model):
    minION = models.ForeignKey(MinION, blank=True, null=True, related_name='minioncontrol')
    owner = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='controlcontrol')
    job = models.CharField(max_length=256, blank=False, null=False)
    custom = models.CharField(max_length=256, blank=True, null=True)
    #setdate = models.DateTimeField(blank=False,null=False)
    #completedate = models.DateTimeField(blank=True,null=True)
    complete = models.BooleanField(default=False)
    created_date = models.DateTimeField(auto_now_add=True)
    modified_date = models.DateTimeField(auto_now=True)  # These last two fields added to enable auto cleanup of event status for a minION incase of disconnect of client.

    class Meta:
        verbose_name = 'MinION Control'
        verbose_name_plural = 'MinION Controls'

    def __str__(self):
        return self.job


class MinIONRun(models.Model):
    run_name = models.CharField(max_length=64)
    run_id = models.CharField(max_length=64)
    is_barcoded = models.BooleanField(default=False)
    owner = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='runs')
    minION = models.ForeignKey(MinION, blank=True, null=True, related_name='minionrun')

    class Meta:
        verbose_name = 'MinION Run'
        verbose_name_plural = 'MinION Runs'

    def __str__(self):
        return self.run_name

    def last_entry(self):
        try:
            return self.RunStats.last().created_date
        except AttributeError:
            return "undefined"

    def last_read(self):
        try:
            #print (self.reads.last().created_date)
            return self.reads.last().created_date
        except AttributeError:
            return "undefined"

    def sample_name(self):
        try:
            return self.RunDetails.last().minKNOW_sample_name
            #return self.RunDetails.minKNOW_sample_name
        except AttributeError:
            return "undefined"

    def minKNOW_flow_cell_id(self):
        try:
            return self.RunDetails.last().minKNOW_flow_cell_id
        except AttributeError:
            return "undefined"

    def minKNOW_version(self):
        try:
            return self.RunDetails.last().minKNOW_version
        except AttributeError:
            return "undefined"

    def barcodes(self):
        barcodes = ['All reads', ]
        for item in FastqRead.objects.filter(run_id=self).values('barcode').distinct():
            barcodes.append(item['barcode'])

        return barcodes


class MinIONStatus(models.Model):
    minION = models.OneToOneField(MinION, related_name='currentdetails')
    minKNOW_status = models.CharField(max_length=64)
    minKNOW_current_script = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_sample_name = models.CharField(max_length=256,blank=True, null=True)
    minKNOW_exp_script_purpose = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_flow_cell_id = models.CharField(max_length=64, blank=True, null=True)
    minKNOW_run_name = models.CharField(max_length=256,blank=True, null=True)
    minKNOW_hash_run_id = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_script_run_id = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_real_sample_rate = models.IntegerField(blank=True,null=True)
    #minKNOW_voltage_offset = models.IntegerField(blank=True, null=True)
    #minKNOW_yield = models.IntegerField(blank=True, null=True)
    minKNOW_asic_id = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_total_drive_space = models.FloatField(blank=True, null=True)
    minKNOW_disk_space_till_shutdown = models.FloatField(blank=True, null=True)
    minKNOW_disk_available = models.FloatField(blank=True,null=True)
    minKNOW_warnings = models.BooleanField(default=False)

    class Meta:
        verbose_name = 'MinION Status'
        verbose_name_plural = 'MinION Status'

    def __str__(self):
        return "{} {}".format(self.minION,self.minKNOW_status)


class MinIONRunStats(models.Model):
    minION = models.ForeignKey(MinION, related_name = 'currentrunstats')
    run_id = models.ForeignKey(MinIONRun, related_name='RunStats')
    sample_time = models.DateTimeField()
    event_yield = models.BigIntegerField(default=0)
    asic_temp = models.FloatField(default=0)
    heat_sink_temp= models.FloatField(default=0)
    voltage_value= models.FloatField(default=0)
    mean_ratio= models.FloatField(default=0)
    open_pore= models.FloatField(default=0)
    in_strand= models.FloatField(default=0)
    multiple = models.IntegerField(default=0)
    unavailable = models.IntegerField(default=0)
    unknown = models.IntegerField(default=0)
    adapter = models.IntegerField(default=0)
    pending_mux_change = models.IntegerField(default=0)
    unclassified = models.IntegerField(default=0)
    below = models.IntegerField(default=0)
    unblocking = models.IntegerField(default=0)
    above = models.IntegerField(default=0)
    good_single = models.IntegerField(default=0)
    saturated = models.IntegerField(default=0)
    inrange = models.IntegerField(default=0)
    strand = models.IntegerField(default=0)
    minKNOW_read_count = models.IntegerField(default=0)
    minKNOW_histogram_values = models.TextField(blank=True, null=True)
    minKNOW_histogram_bin_width = models.IntegerField(default=900)
    created_date = models.DateTimeField(auto_now_add=True)

    class Meta:
        verbose_name = 'MinION Run Stats'
        verbose_name_plural = 'MinION Run Stats'

    def __str__(self):
        return "{} {} {}".format(self.minION,self.run_id,self.sample_time)

    def occupancy(self):
        if (self.strand > 0 and self.inrange > 0):
            occupancy = (self.strand / (self.strand + self.inrange)) * 100
        else:
            occupancy = 0
        return occupancy


class MinIONRunStatus(models.Model):
    minION = models.ForeignKey(MinION, related_name='currentrundetails')
    #minKNOW_status = models.CharField(max_length=64)
    minKNOW_current_script = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_sample_name = models.CharField(max_length=256,blank=True, null=True)
    minKNOW_exp_script_purpose = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_flow_cell_id = models.CharField(max_length=64, blank=True, null=True)
    minKNOW_version = models.CharField(max_length=64, blank=True, null=True)
    minKNOW_run_name = models.CharField(max_length=256,blank=True, null=True)
    run_id = models.ForeignKey(MinIONRun , related_name='RunDetails')
    minKNOW_hash_run_id = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_script_run_id = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_real_sample_rate = models.IntegerField(blank=True,null=True)
    #minKNOW_voltage_offset = models.IntegerField(blank=True, null=True)
    #minKNOW_yield = models.IntegerField(blank=True, null=True)
    minKNOW_asic_id = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_start_time = models.DateTimeField(blank=True,null=True)
    minKNOW_colours_string = models.TextField(blank=True,null=True)
    #minKNOW_total_drive_space = models.FloatField(blank=True, null=True)
    #minKNOW_disk_space_till_shutdown = models.FloatField(blank=True, null=True)
    #minKNOW_warnings = models.BooleanField(default=False)

    class Meta:
        verbose_name = 'MinION Run Status'
        verbose_name_plural = 'MinION Run Status'

    def __str__(self):
        return "{} {} {}".format(self.minION,self.minKNOW_current_script, self.run_id)

    def minION_name(self):
        return self.minION.minION_name


class MinIONEventType(models.Model):
    name = models.CharField(max_length=64)

    class Meta:
        verbose_name = 'MinION Event Type'
        verbose_name_plural = 'MinION Event Types'

    def __str__(self):
        return self.name


class MinIONEvent(models.Model):
    computer_name = models.CharField(max_length=256)
    minION = models.ForeignKey(MinION, related_name='events')
    event = models.ForeignKey(MinIONEventType)
    datetime =models.DateTimeField()
    created_date = models.DateTimeField(auto_now_add=True)
    modified_date = models.DateTimeField(auto_now=True)  #These last two fields added to enable auto cleanup of event status for a minION incase of disconnect of client.

    class Meta:
        verbose_name = 'MinION Event'
        verbose_name_plural = 'MinION Events'

    def __str__(self):
        return "{} {} {} {}".format(self.computer_name, self.minION, self.event, self.datetime)


class MinIONScripts(models.Model):
    minION = models.ForeignKey(MinION, related_name='scripts')
    identifier = models.CharField(max_length=256)
    name = models.CharField(max_length=256)
    experiment_type = models.CharField(max_length=256, blank=True, null=True)
    base_calling = models.NullBooleanField(blank=True, null=True)
    flow_cell = models.CharField(max_length=256, blank=True, null=True)
    kit = models.CharField(max_length=256, blank=True, null=True)

    class Meta:
        verbose_name = 'MinION Script'
        verbose_name_plural = 'MinION Scripts'

    def __str__(self):
        return "{} {} {}".format(self.minION, self.name, self.identifier)


class Read(models.Model):
    pass


class FastqReadType(models.Model):
    name = models.CharField(max_length=16)

    class Meta:
        verbose_name = 'FASTQ Read Type'
        verbose_name_plural = 'FASTQ Read Types'

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
    created_date = models.DateTimeField(auto_now_add=True)
    modified_date = models.DateTimeField(auto_now=True)

    class Meta:
        verbose_name = 'FASTQ Read'
        verbose_name_plural = 'FASTQ Read'

    def __str__(self):
        return self.read_id


class RunSummary(models.Model):
    run_id = models.ForeignKey(MinIONRun,on_delete=models.CASCADE, related_name='runsummaries')
    total_length = models.BigIntegerField(default=0)
    read_count = models.IntegerField(default=0)
    type = models.ForeignKey(FastqReadType)
    max_length = models.IntegerField(default=0)
    min_length = models.IntegerField(default=0)
    pass_length = models.BigIntegerField(default=0)
    pass_max_length = models.IntegerField(default=0)
    pass_min_length = models.IntegerField(default=0)
    pass_count = models.IntegerField(default=0)

    class Meta:
        verbose_name = 'Run Summary'
        verbose_name_plural = 'Run Summaries'

    def __str__(self):
        return "{} {} {} {}".format(self.run_id, self.total_length, self.read_count, self.type)


class RunSummaryBarCode(models.Model):
    run_id = models.ForeignKey(MinIONRun,on_delete=models.CASCADE, related_name='runsummariesbarcodes')
    total_length = models.BigIntegerField(default=0)
    read_count = models.IntegerField(default=0)
    type = models.ForeignKey(FastqReadType)
    barcode = models.CharField(max_length=32)
    pass_length = models.BigIntegerField(default=0)
    pass_max_length = models.IntegerField(default=0)
    pass_min_length = models.IntegerField(default=0)
    pass_count = models.IntegerField(default=0)
    max_length = models.IntegerField(default=0)
    min_length = models.IntegerField(default=0)

    class Meta:
        verbose_name = 'Run Summary Barcode'
        verbose_name_plural = 'Run Summary Barcodes'

    def __str__(self):
        return "{} {} {} {} {}".format(self.run_id, self.total_length, self.read_count, self.type, self.barcode)

class ChannelSummary(models.Model):
    run_id = models.ForeignKey(MinIONRun, on_delete=models.CASCADE)
    channel_number = models.IntegerField()
    read_count = models.BigIntegerField()
    read_length = models.BigIntegerField()

    def __str__(self):
        return "{} {} {}".format(self.run_id,self.channel_number,self.read_count)

class RunStatistic(models.Model):
    run_id = models.ForeignKey(MinIONRun, on_delete=models.CASCADE)
    sample_time = models.DateTimeField()
    total_length = models.BigIntegerField(default=0)
    max_length = models.IntegerField(default=0)
    min_length = models.IntegerField(default=0)
    read_count = models.IntegerField(default=0)
    pass_length = models.BigIntegerField(default=0)
    pass_max_length = models.IntegerField(default=0)
    pass_min_length = models.IntegerField(default=0)
    pass_count = models.IntegerField(default=0)
    type = models.ForeignKey(FastqReadType)

    class Meta:
        verbose_name = 'Run Statistics'
        verbose_name_plural = 'Run Statistics'

    def __str__(self):
        return "{} {} {}".format(self.run_id, self.sample_time, self.type)


class RunStatisticBarcode(models.Model):
    run_id = models.ForeignKey(MinIONRun, on_delete=models.CASCADE)
    sample_time = models.DateTimeField()
    total_length = models.BigIntegerField(default=0)
    read_count = models.IntegerField(default=0)
    max_length = models.IntegerField(default=0)
    min_length = models.IntegerField(default=0)
    pass_length = models.BigIntegerField(default=0)
    pass_max_length = models.IntegerField(default=0)
    pass_min_length = models.IntegerField(default=0)
    pass_count = models.IntegerField(default=0)
    type = models.ForeignKey(FastqReadType)
    barcode = models.CharField(max_length=32)

    class Meta:
        verbose_name = 'Run Statistics Barcode'
        verbose_name_plural = 'Run Statistics Barcodes'

    def __str__(self):
        return "{} {} {} {}".format(self.run_id, self.sample_time, self.type, self.barcode)


class MinIONmessages(models.Model):
    minION = models.ForeignKey(MinION, related_name='messages')
    run_id = models.ForeignKey(MinIONRun, related_name='runmessages', blank=True, null=True)
    minKNOW_message = models.CharField(max_length=256)
    minKNOW_identifier = models.CharField(max_length=256)
    minKNOW_severity = models.CharField(max_length=64)
    minKNOW_message_timestamp = models.DateTimeField()

    class Meta:
        unique_together = (("minION", "run_id", "minKNOW_message_timestamp"))
        verbose_name = 'MinION Message'
        verbose_name_plural = 'MinION Messages'

    def __str__(self):
        return "{} {} {} {}".format(self.minION, self.minKNOW_message, self.minKNOW_severity, self.minKNOW_message_timestamp)


##### Job Management System

class Job(models.Model):
    jobname = models.CharField(max_length=256)

    def __str__(self):
        return "{}".format(self.jobname)

class JobMaster(models.Model):
    run_id = models.ForeignKey(MinIONRun, related_name='runjobs')
    job_name = models.ForeignKey(Job, related_name='taskname')
    var1 = models.CharField(max_length=256,blank=True, null=True)
    var2 = models.CharField(max_length=256,blank=True, null=True)
    var3 = models.CharField(max_length=256,blank=True, null=True)
    complete = models.BooleanField(default=False)
    running = models.BooleanField(default=False)


    def __str__(self):
        return "{} {}".format(self.run_id,self.job_name)

##### Alignment Management


##### This class replicates original minoTour design but perhaps it is better just to store the sam file as a single row? Lets try that
"""
class SamStore(models.Model):
    run_id = models.ForeignKey(MinIONRun, related_name='runalignments')
    read_id = models.ForeignKey(FastqRead, related_name='readalignment')
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

class SamStore(models.Model):
    run_id = models.ForeignKey(MinIONRun, related_name='runalignments')
    read_id = models.ForeignKey(FastqRead, related_name='readalignment')
    samline = models.TextField()

    def __str__(self):
        return "{} {}".format(self.run_id,self.read_id)

@receiver(post_save, sender=settings.AUTH_USER_MODEL)
def create_auth_token(sender, instance=None, created=False, **kwargs):
    if created:
        Token.objects.create(user=instance)


@receiver(post_save, sender=FastqRead)
def update_global_state(instance, sender, **kwargs):
    ipn_obj = instance

    barcode = ipn_obj.barcode

    tm = ipn_obj.start_time

    tm = tm - datetime.timedelta(minutes=(tm.minute % 1) - 1,
                                 seconds=tm.second,
                                 microseconds=tm.microsecond)

    obj1,created1 = RunSummaryBarCode.objects.update_or_create(
        run_id=ipn_obj.run_id, type=ipn_obj.type, barcode='All reads'
    )
    update_sum_stats(obj1, ipn_obj)

    # all reads and barcodes are saved on RunStatisticBarcode
    obj3, created3 = RunStatisticBarcode.objects.update_or_create(
        run_id=ipn_obj.run_id, type=ipn_obj.type, barcode='All reads', sample_time=tm
    )
    update_sum_stats(obj3, ipn_obj)

    if barcode is not None and barcode != '':
        obj2,created2 = RunSummaryBarCode.objects.update_or_create(
            run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=barcode
        )
        update_sum_stats(obj2, ipn_obj)

        obj3, created3 = RunStatisticBarcode.objects.update_or_create(
            run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=barcode, sample_time=tm
        )
        update_sum_stats(obj3, ipn_obj)


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
