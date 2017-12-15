import datetime

import pytz
from django.conf import settings
from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver
from rest_framework.authtoken.models import Token

from reference.models import ReferenceInfo


class FlowCell(models.Model):
    name = models.CharField(
        max_length=256,
        blank=True,
        null=True,
    )
    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        related_name='flowcells'
    )

    size = models.IntegerField(
        blank=True,
        null=True,
    )

    def __str__(self):
        return "{} {}".format(self.name, self.id)

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


class UserOptions(models.Model):
    owner = models.OneToOneField(settings.AUTH_USER_MODEL, related_name='extendedopts')
    twitterhandle = models.CharField(max_length=64)
    tweet = models.BooleanField(default=False)
    email = models.BooleanField(default=False)

    def __str__(self):
        return "{}".format(str(self.owner))


class MinIONRun(models.Model):
    run_name = models.CharField(
        max_length=256
    )

    run_id = models.CharField(
        max_length=64
    )

    is_barcoded = models.BooleanField(
        default=False
    )

    has_fastq = models.BooleanField(
        default=True
    )

    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        related_name='runs'
    )

    minION = models.ForeignKey(
        MinION,
        blank=True,
        null=True,
        related_name='minionrun'
    )

    active = models.BooleanField(
        default=False
    )

    to_delete = models.BooleanField(
        default=False
    )

    class Meta:
        verbose_name = 'MinION Run'
        verbose_name_plural = 'MinION Runs'

    def __str__(self):
        return "{} - {}".format(self.run_name, self.run_id)

    def last_entry(self):
        try:
            return self.RunStats.last().created_date
        except AttributeError:
            #return "undefined"
            olddate = datetime.datetime(1, 1, 1, 1, 1, 1, 1, pytz.utc)
            return olddate


    def start_time(self):
        try:
            return self.RunDetails.last().minKNOW_start_time
        except AttributeError:
            return "undefined"

    def last_read(self):
        try:
            return self.reads.last().created_date
        except AttributeError:
            #return "undefined"
            olddate = datetime.datetime(1, 1, 1, 1, 1, 1, 1, pytz.utc)
            return olddate

    def sample_name(self):
        try:
            return self.RunDetails.last().minKNOW_sample_name
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

    def max_channel(self):
        try:
            return self.runchansum.order_by('channel_number').last().channel_number
        except AttributeError:
            return "undefined"

    def flowcell_type(self):
        try:
            max_channel = self.max_channel()
            if max_channel != 'undefined':
                if int(max_channel) <= 128:
                    return 128
                elif 128 < int(max_channel) <= 512:
                    return 512
                elif 512 < int(max_channel) <=3000:
                    return 3000
            return 512
        except AttributeError:
            return "undefined"

    #def barcodes(self):
    #    barcodes = []
    #    for item in Barcode.objects.filter(run=self):
    #        barcodes.append(item.name)
    #
    #    return barcodes


class BarcodeGroup(models.Model):
    flowcell = models.ForeignKey(
        FlowCell,
        related_name='flowcellbarcode'
    )

    name = models.CharField(
        max_length=32
    )

    def __str__(self):
        return "{} {}".format(self.flowcell,self.name)


class Barcode(models.Model):
    run = models.ForeignKey(
        MinIONRun,
        on_delete=models.CASCADE,
        related_name = 'barcodes'
    )

    barcodegroup = models.ForeignKey(
        BarcodeGroup,
        on_delete=models.CASCADE,
        null=True,blank=True,
        related_name = 'barcodegroup'
    )

    name = models.CharField(
        max_length=32
    )

    def __str__(self):
        return "{} {} {}".format(self.run, self.run.run_id, self.name)




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


    ## This is something to look at for optimisation
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
    minKNOW_computer = models.TextField(max_length=128,blank=True,null=True)
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


class FastqReadType(models.Model):
    name = models.CharField(
        max_length=16
    )

    class Meta:
        verbose_name = 'FASTQ Read Type'
        verbose_name_plural = 'FASTQ Read Types'

    def __str__(self):
        return self.name


class FastqRead(models.Model):
    run_id = models.ForeignKey(
        MinIONRun,
        on_delete=models.CASCADE,
        related_name='reads'
    )

    read_id = models.CharField(
        max_length=96
    )

    read = models.IntegerField(

    )

    channel = models.IntegerField(

    )

    barcode = models.ForeignKey(
        Barcode,
        #on_delete=models.CASCADE,
        related_name='reads',
        null=True
    )

    sequence_length = models.BigIntegerField(
        null=True,
        blank=True
    )

    quality_average = models.DecimalField(
        decimal_places=2,
        max_digits=5,
        null=True,
        blank=True
    )

    is_pass = models.BooleanField(

    ) # pass = true, fail = false

    type = models.ForeignKey(
        FastqReadType
    )

    start_time = models.DateTimeField(

    )

    created_date = models.DateTimeField(
        auto_now_add=True
    )

    modified_date = models.DateTimeField(
        auto_now=True
    )

    class Meta:
        verbose_name = 'FASTQ Read'
        verbose_name_plural = 'FASTQ Read'

    def __str__(self):
        return self.read_id



class FastqReadExtra(models.Model):
    fastqread = models.OneToOneField(
        FastqRead,
        on_delete=models.CASCADE,
        related_name='extra'
    )

    sequence = models.TextField(
        blank=True,
        null=True
    )

    quality = models.TextField(
        blank=True,
        null=True,
    )

    def __str__(self):
        return self.fastqread

class RunSummary(models.Model):
    run_id = models.ForeignKey(MinIONRun,on_delete=models.CASCADE, related_name='runsummaries')
    total_length = models.BigIntegerField(default=0)
    read_count = models.IntegerField(default=0)
    type = models.ForeignKey(FastqReadType)
    quality_sum = models.DecimalField(
        decimal_places=2,
        max_digits=12,
        default=0
    )
    max_length = models.IntegerField(default=0)
    min_length = models.IntegerField(default=0)
    pass_quality_sum = models.DecimalField(
        decimal_places=2,
        max_digits=12,
        default=0
    )
    pass_length = models.BigIntegerField(default=0)
    pass_max_length = models.IntegerField(default=0)
    pass_min_length = models.IntegerField(default=0)
    pass_count = models.IntegerField(default=0)

    class Meta:
        verbose_name = 'Run Summary'
        verbose_name_plural = 'Run Summaries'

    def __str__(self):
        return "{} {} {} {}".format(self.run_id, self.total_length, self.read_count, self.type)


class RunSummaryBarcode(models.Model):
    run_id = models.ForeignKey(
        MinIONRun,
        on_delete=models.CASCADE,
        related_name='runsummariesbarcodes'
    )

    total_length = models.BigIntegerField(
        default=0
    )

    quality_sum = models.DecimalField(
        decimal_places=2,
        max_digits=12,
        default=0
    )

    read_count = models.IntegerField(
        default=0
    )

    type = models.ForeignKey(
        FastqReadType
    )

    barcode = models.ForeignKey(
        Barcode,
        on_delete=models.CASCADE,
        related_name='runsummaries',
        null=True
    )

    pass_quality_sum = models.DecimalField(
        decimal_places=2,
        max_digits=12,
        default=0
    )

    pass_length = models.BigIntegerField(
        default=0
    )

    pass_max_length = models.IntegerField(
        default=0
    )

    pass_min_length = models.IntegerField(
        default=0
    )

    pass_count = models.IntegerField(
        default=0
    )

    max_length = models.IntegerField(
        default=0
    )

    min_length = models.IntegerField(
        default=0
    )

    channel_presence = models.CharField(
        max_length=3000,
        default='0' * 3000
    )

    class Meta:
        verbose_name = 'Run Summary Barcode'
        verbose_name_plural = 'Run Summary Barcodes'

    def __str__(self):
        return "{} {} {} {} {}".format(
            self.run_id,
            self.total_length,
            self.read_count,
            self.type,
            self.barcode
        )

    def number_active_channels(self):
        return len(self.channel_presence.replace('0', ''))


class HistogramSummary(models.Model):
    run_id = models.ForeignKey(MinIONRun, on_delete=models.CASCADE)
    barcode = models.ForeignKey(Barcode, on_delete=models.CASCADE)
    read_type = models.ForeignKey(FastqReadType)
    bin_width = models.BigIntegerField()
    read_count = models.BigIntegerField(default=0)
    read_length = models.BigIntegerField(default=0)

    def __str__(self):
        return "{} {} {}".format(self.run_id,self.read_type,self.bin_width)


class ChannelSummary(models.Model):
    run_id = models.ForeignKey(MinIONRun, on_delete=models.CASCADE, related_name='runchansum')
    channel_number = models.IntegerField()
    read_count = models.BigIntegerField(default=0)
    read_length = models.BigIntegerField(default=0)

    def __str__(self):
        return "{} {} {}".format(self.run_id,self.channel_number,self.read_count)


class RunStatistic(models.Model):
    run_id = models.ForeignKey(MinIONRun, on_delete=models.CASCADE)
    sample_time = models.DateTimeField()
    total_length = models.BigIntegerField(default=0)
    quality_sum = models.DecimalField(
        decimal_places=2,
        max_digits=12,
        default=0
    )
    max_length = models.IntegerField(default=0)
    min_length = models.IntegerField(default=0)
    read_count = models.IntegerField(default=0)
    pass_quality_sum = models.DecimalField(
        decimal_places=2,
        max_digits=12,
        default=0
    )
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
    run_id = models.ForeignKey(
        MinIONRun,
        on_delete=models.CASCADE,
        related_name='runstatbarc',
    )

    sample_time = models.DateTimeField(

    )

    total_length = models.BigIntegerField(
        default=0
    )

    read_count = models.IntegerField(
        default=0
    )

    max_length = models.IntegerField(
        default=0
    )

    min_length = models.IntegerField(
        default=0
    )

    pass_length = models.BigIntegerField(
        default=0
    )

    pass_max_length = models.IntegerField(
        default=0
    )

    pass_min_length = models.IntegerField(
        default=0
    )

    pass_count = models.IntegerField(
        default=0
    )
    quality_sum = models.DecimalField(
        decimal_places=2,
        max_digits=12,
        default=0
    )
    pass_quality_sum = models.DecimalField(
        decimal_places=2,
        max_digits=12,
        default=0
    )

    type = models.ForeignKey(
        FastqReadType
    )

    barcode = models.ForeignKey(
        Barcode,
        on_delete=models.CASCADE,
        related_name='runstatistics',
        null=True
    )

    channel_presence = models.CharField(
        max_length=3000,
        default='0' * 3000
    )

    class Meta:
        verbose_name = 'Run Statistics Barcode'
        verbose_name_plural = 'Run Statistics Barcodes'

    def __str__(self):
        return "{} {} {} {}".format(
            self.run_id,
            self.sample_time,
            self.type,
            self.barcode
        )

    def number_active_channels(self):
        return len(self.channel_presence.replace('0', ''))


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
        return "{} {} {} {}".format(
            self.minION, self.minKNOW_message, self.minKNOW_severity, self.minKNOW_message_timestamp)




#@receiver(post_save, sender=MinIONRun)
def create_all_reads_barcode(sender, instance=None, created=False, **kwargs):
    if created:
        #barcodegroup, created = BarcodeGroup.objects.get_or_create(flowcell=instance.flowcellrun.first().flowcell,
        #                                                           name='All reads')

        Barcode.objects.get_or_create(run=instance, name='All reads')


        #barcodegroup, created = BarcodeGroup.objects.get_or_create(flowcell=instance.flowcellrun.first().flowcell,
        #                                                           name='No barcode')

        Barcode.objects.get_or_create(run=instance, name='No barcode')




@receiver(post_save, sender=settings.AUTH_USER_MODEL)
def create_auth_token(sender, instance=None, created=False, **kwargs):
    if created:
        Token.objects.create(user=instance)


#@receiver(post_save, sender=FastqRead)
def update_global_state(instance): ##, sender, **kwargs):
#    print ("trying to run")
#    update_global_state_task.delay(instance)

#@task()
#def update_global_state_task(instance):
#    print ("update task running")
    ipn_obj = instance

    barcode = ipn_obj.barcode

    barcode_all_reads = ipn_obj.run_id.barcodes.filter(name='All reads').first()

    tm = ipn_obj.start_time

    tm = tm - datetime.timedelta(minutes=(tm.minute % 1) - 1,
                                 seconds=tm.second,
                                 microseconds=tm.microsecond)

    obj1, created1 = RunSummaryBarcode.objects.update_or_create(
        run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=barcode_all_reads
    )
    update_sum_stats(obj1, ipn_obj)

    # all reads and barcodes are saved on RunStatisticBarcode
    obj3, created3 = RunStatisticBarcode.objects.update_or_create(
        run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=barcode_all_reads, sample_time=tm
    )
    update_sum_stats(obj3, ipn_obj)

    if barcode is not None and barcode.name != '':
        obj2,created2 = RunSummaryBarcode.objects.update_or_create(
            run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=barcode
        )
        update_sum_stats(obj2, ipn_obj)

        obj3, created3 = RunStatisticBarcode.objects.update_or_create(
            run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=barcode, sample_time=tm
        )
        update_sum_stats(obj3, ipn_obj)


def update_sum_stats(obj, ipn_obj):

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

    #
    # channel_presence is a 512 characters length string containing 0
    # for each channel (id between 1 - 512) seen, we set the position on
    # the string to 1. afterwards, we can calculate the number of active
    # channels in a particular minute.
    #
    channel = ipn_obj.channel
    channel_sequence = obj.channel_presence
    channel_sequence_list = list(channel_sequence)
    channel_sequence_list[channel] = '1'
    obj.channel_presence = ''.join(channel_sequence_list)

    obj.read_count += 1
    obj.save()



class GroupRunType(models.Model):

    type_name = models.CharField(
        max_length=256
    )

    def __str__(self):
        return "{}".format(self.type_name)


class JobType(models.Model):

    name = models.CharField(
        max_length=256
    )

    description = models.TextField(
        max_length=256,
        blank = True,
        null = True
    )

    long_description = models.TextField(
        blank=True,
        null=True
    )

    reference = models.BooleanField(
        default=False,
    )

    transcriptome = models.BooleanField(
        default=False,
    )

    readcount = models.BooleanField(
        default=False,
    )
    private = models.BooleanField(
        default=True
    )

    type_name = models.ManyToManyField(GroupRunType)


    def __str__(self):
        return "{}".format(self.name)




class JobMaster(models.Model):

    run = models.ForeignKey(
        MinIONRun,
        on_delete=models.CASCADE,
        related_name='runjobs',
        null=True,
        blank=True
    )

    flowcell = models.ForeignKey(
        FlowCell,
        on_delete=models.CASCADE,
        related_name='flowcelljobs',
        null=True,
        blank=True
    )

    job_type = models.ForeignKey(
        JobType,
        related_name='taskname'
    )

    reference = models.ForeignKey(
        ReferenceInfo,
        related_name='referencejob',
        null=True,
        blank=True
    )

    last_read = models.CharField(
        max_length=256,
        blank=True,
        null=True
    )

    tempfile_name = models.CharField(
        max_length=256,
        blank=True,
        null=True
    )

    read_count = models.BigIntegerField(
        default=0
    )

    complete = models.BooleanField(
        default=False
    )

    running = models.BooleanField(
        default=False
    )

    def __str__(self):
        if self.run is not None:
            return "{} {} {}".format(self.run, self.job_type, self.run.id)
        else:
            return "{} {} {}".format(self.flowcell, self.job_type, self.flowcell.id)






class FlowCellRun(models.Model):
    flowcell = models.ForeignKey(
        FlowCell,
        related_name='flowcelldetails'
    )

    run = models.ForeignKey(
        MinIONRun,
        on_delete=models.CASCADE,
        related_name='flowcellrun'
    )

    def name(self):
        return self.flowcell.name

    def run_name(self):
        return self.run.run_name

    def barcodes(self):
        return self.run.barcodes

    def __str__(self):
        return "{} {}".format(self.flowcell, self.run)

    def last_entry(self):
        try:
            return self.run.runstatbarc.last().sample_time
            #return self.RunStats.last().created_date
        except AttributeError:
            #return "undefined"
            olddate = datetime.datetime(1, 1, 1, 1, 1, 1, 1, pytz.utc)
            return olddate

    def start_time(self):
        try:
            return self.run.runstatbarc.first().sample_time
            #return self.RunDetails.last().minKNOW_start_time
        except AttributeError:
            return "undefined"

    def last_read(self):
        try:
            return self.run.runstatbarc.last().sample_time
            #return self.reads.last().created_date
        except AttributeError:
            #return "undefined"
            olddate = datetime.datetime(1, 1, 1, 1, 1, 1, 1, pytz.utc)
            return olddate

    def sample_name(self):
        try:
            return self.run.minKNOW_sample_name()
        except AttributeError:
            return "undefined"

    def minKNOW_flow_cell_id(self):
        try:
            return self.run.minKNOW_flow_cell_id()
        except AttributeError:
            return "undefined"

    def minKNOW_version(self):
        try:
            return self.run.minKNOW_version()
        except AttributeError:
            return "undefined"

    def active(self):
        return self.run.active



@receiver(post_save, sender=FlowCellRun)
def create_all_reads_barcode_flowcellrun(sender, instance=None, created=False, **kwargs):
    if created:
        barcodegroup, created = BarcodeGroup.objects.get_or_create(flowcell=instance.flowcell,
                                                                   name='All reads')

        Barcode.objects.update_or_create(run=instance.run, name='All reads', barcodegroup=barcodegroup)

        barcodegroup, created = BarcodeGroup.objects.get_or_create(flowcell=instance.flowcell,
                                                                   name='No barcode')

        Barcode.objects.update_or_create(run=instance.run, name='No barcode', barcodegroup=barcodegroup)