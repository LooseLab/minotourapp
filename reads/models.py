import datetime
import hashlib

import pytz
from django.conf import settings
from django.contrib.auth.models import User
from django.db import models
from django.db.models.signals import post_save, m2m_changed
from django.dispatch import receiver
from rest_framework.authtoken.models import Token

from devices.models import Flowcell, MinION


def create_flowcell(owner_id, name=""):

    if name == "":

        text = '{}'.format(datetime.datetime.now())
        h = hashlib.md5()
        h.update(text.encode('utf-8'))
        name = "flowcell-{}".format(h.hexdigest())[0:15]

    owner = User.objects.get(pk=owner_id)

    return Flowcell.objects.create(name=name, owner=owner)


class GroupRun(models.Model):

    MINION = 'MINION'
    GRIDION = 'GRIDION'
    PROMETHION = 'PROMETHION'

    NANOPORE_DEVICES = (

        (MINION, 'MinION'),
        (GRIDION, 'GridION'),
        (PROMETHION, 'PromethION')
    )

    owner = models.ForeignKey(

        settings.AUTH_USER_MODEL,
        related_name='groupruns'
    )

    name = models.CharField(

        max_length = 32
    )

    device = models.CharField(

        max_length = 10,
        choices = NANOPORE_DEVICES,
        default = MINION
    )

    def __str__(self):
        return self.name

    # def barcodes(self):
    #
    #     barcode_list = []
    #
    #     run_list = Run.objects.filter(groupruns=self)
    #
    #     for run in run_list:
    #
    #         for barcode in run.barcodes.all():
    #
    #             barcode_list.append(barcode)
    #
    #     return barcode_list


class MinIONControl(models.Model):
    minION = models.ForeignKey(MinION, blank=True, null=True, related_name='minioncontrol')
    owner = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='controlcontrol')
    job = models.CharField(max_length=256, blank=False, null=False)
    custom = models.CharField(max_length=256, blank=True, null=True)
    # setdate = models.DateTimeField(blank=False,null=False)
    # completedate = models.DateTimeField(blank=True,null=True)
    complete = models.BooleanField(default=False)
    created_date = models.DateTimeField(auto_now_add=True)
    modified_date = models.DateTimeField(
        auto_now=True)  # These last two fields added to enable auto cleanup of event status for a minION incase of disconnect of client.

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


class Run(models.Model):

    minion = models.ForeignKey(

        MinION,
        blank=True,
        null=True,
        related_name='runs'
    )

    flowcell = models.ForeignKey(

        Flowcell,
        related_name='runs',
        null=True,
        blank=True,
    )

    owner = models.ForeignKey(

        settings.AUTH_USER_MODEL,
        related_name='runs'
    )

    groupruns = models.ManyToManyField(

        GroupRun,
        blank = True,
        related_name='runs',
    )

    name = models.CharField(

        max_length=256
    )

    runid = models.CharField(

        max_length=64
    )

    is_barcoded = models.BooleanField(

        default=False
    )

    has_fastq = models.BooleanField(

        default=True
    )

    active = models.BooleanField(

        default=False
    )

    to_delete = models.BooleanField(

        default=False
    )

    class Meta:

        verbose_name = 'Run'
        verbose_name_plural = 'Runs'

    """def save(self, *args, **kwargs):

        if self.flowcell is None:
            self.flowcell = create_flowcell(self.owner.id)

        super(Run, self).save(*args, **kwargs)"""


    def __str__(self):

        return "{} - {}".format(self.name, self.runid)

    def flowcell_name(self):

        try:
            return self.flowcellrun.last().flowcell.name
        except AttributeError:
            return "undefined"

    def last_entry(self):

        try:
            return self.RunStats.last().created_date
        except AttributeError:
            # return "undefined"
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
            # return "undefined"
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
            return self.runchansum.order_by('channel').last().channel
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
                elif 512 < int(max_channel) <= 3000:
                    return 3000
            return 512
        except AttributeError:
            return "undefined"


class BarcodeGroup(models.Model):

    flowcell = models.ForeignKey(
        Flowcell,
        related_name='flowcellbarcode'
    )

    name = models.CharField(
        max_length=32
    )

    def __str__(self):
        return "{} {}".format(self.flowcell, self.name)


class GroupBarcode(models.Model):

    grouprun = models.ForeignKey(

        GroupRun,
        related_name='groupbarcodes'
    )

    name = models.CharField(

        max_length=32
    )

    def __str__(self):
        return "{} {}".format(self.grouprun, self.name)


class Barcode(models.Model):

    run = models.ForeignKey(

        Run,
        on_delete=models.CASCADE,
        related_name='barcodes'
    )

    barcodegroup = models.ForeignKey(

        BarcodeGroup,
        on_delete=models.CASCADE,
        null=True, blank=True,
        related_name='barcodegroup'
    )

    groupbarcodes = models.ManyToManyField(

        GroupBarcode,

        blank=True,
        related_name='barcodes'
    )

    name = models.CharField(

        max_length=32
    )

    def __str__(self):

        return "{} {} {}".format(self.run, self.run.runid, self.name)


class MinIONStatus(models.Model):
    minION = models.OneToOneField(MinION, related_name='currentdetails')
    minKNOW_status = models.CharField(max_length=64)
    minKNOW_current_script = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_sample_name = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_exp_script_purpose = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_flow_cell_id = models.CharField(max_length=64, blank=True, null=True)
    minKNOW_run_name = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_hash_run_id = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_script_run_id = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_real_sample_rate = models.IntegerField(blank=True, null=True)
    # minKNOW_voltage_offset = models.IntegerField(blank=True, null=True)
    # minKNOW_yield = models.IntegerField(blank=True, null=True)
    minKNOW_asic_id = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_total_drive_space = models.FloatField(blank=True, null=True)
    minKNOW_disk_space_till_shutdown = models.FloatField(blank=True, null=True)
    minKNOW_disk_available = models.FloatField(blank=True, null=True)
    minKNOW_warnings = models.BooleanField(default=False)

    class Meta:
        verbose_name = 'MinION Status'
        verbose_name_plural = 'MinION Status'

    def __str__(self):
        return "{} {}".format(self.minION, self.minKNOW_status)


class MinIONRunStats(models.Model):
    minION = models.ForeignKey(MinION, related_name='currentrunstats')
    run_id = models.ForeignKey(Run, related_name='RunStats')
    sample_time = models.DateTimeField()
    event_yield = models.BigIntegerField(default=0)
    asic_temp = models.FloatField(default=0)
    heat_sink_temp = models.FloatField(default=0)
    voltage_value = models.FloatField(default=0)
    mean_ratio = models.FloatField(default=0)
    open_pore = models.FloatField(default=0)
    in_strand = models.FloatField(default=0)
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
        return "{} {} {}".format(self.minION, self.run_id, self.sample_time)

    ## This is something to look at for optimisation
    def occupancy(self):
        if (self.strand > 0 and self.inrange > 0):
            occupancy = (self.strand / (self.strand + self.inrange)) * 100
        else:
            occupancy = 0
        return occupancy


class MinIONRunStatus(models.Model):
    minION = models.ForeignKey(MinION, related_name='currentrundetails')
    # minKNOW_status = models.CharField(max_length=64)
    minKNOW_current_script = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_sample_name = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_exp_script_purpose = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_flow_cell_id = models.CharField(max_length=64, blank=True, null=True)
    minKNOW_version = models.CharField(max_length=64, blank=True, null=True)
    minKNOW_run_name = models.CharField(max_length=256, blank=True, null=True)
    run_id = models.ForeignKey(Run, related_name='RunDetails')
    minKNOW_hash_run_id = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_script_run_id = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_real_sample_rate = models.IntegerField(blank=True, null=True)
    # minKNOW_voltage_offset = models.IntegerField(blank=True, null=True)
    # minKNOW_yield = models.IntegerField(blank=True, null=True)
    minKNOW_asic_id = models.CharField(max_length=256, blank=True, null=True)
    minKNOW_start_time = models.DateTimeField(blank=True, null=True)
    minKNOW_colours_string = models.TextField(blank=True, null=True)
    minKNOW_computer = models.TextField(max_length=128, blank=True, null=True)

    # minKNOW_total_drive_space = models.FloatField(blank=True, null=True)
    # minKNOW_disk_space_till_shutdown = models.FloatField(blank=True, null=True)
    # minKNOW_warnings = models.BooleanField(default=False)

    class Meta:
        verbose_name = 'MinION Run Status'
        verbose_name_plural = 'MinION Run Status'

    def __str__(self):
        return "{} {} {}".format(self.minION, self.minKNOW_current_script, self.run_id)

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
    datetime = models.DateTimeField()
    created_date = models.DateTimeField(auto_now_add=True)
    modified_date = models.DateTimeField(
        auto_now=True)  # These last two fields added to enable auto cleanup of event status for a minION incase of disconnect of client.

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

    run = models.ForeignKey(
        Run,
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
        # on_delete=models.CASCADE,
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

    )  # pass = true, fail = false

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
        related_name='fastqreadextra'
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


class MinIONmessages(models.Model):
    minION = models.ForeignKey(MinION, related_name='messages')
    run_id = models.ForeignKey(Run, related_name='runmessages', blank=True, null=True)
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


# @receiver(post_save, sender=MinIONRun)
def create_all_reads_barcode(sender, instance=None, created=False, **kwargs):
    if created:
        # barcodegroup, created = BarcodeGroup.objects.get_or_create(flowcell=instance.flowcellrun.first().flowcell,
        #                                                           name='All reads')

        Barcode.objects.get_or_create(run=instance, name='All reads')

        # barcodegroup, created = BarcodeGroup.objects.get_or_create(flowcell=instance.flowcellrun.first().flowcell,
        #                                                           name='No barcode')

        Barcode.objects.get_or_create(run=instance, name='No barcode')


@receiver(post_save, sender=settings.AUTH_USER_MODEL)
def create_auth_token(sender, instance=None, created=False, **kwargs):
    if created:
        Token.objects.create(user=instance)


class FlowCellRun(models.Model):
    flowcell = models.ForeignKey(
        Flowcell,
        related_name='flowcelldetails'
    )

    run = models.ForeignKey(
        Run,
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
            # return self.RunStats.last().created_date
        except AttributeError:
            # return "undefined"
            olddate = datetime.datetime(1, 1, 1, 1, 1, 1, 1, pytz.utc)
            return olddate

    def start_time(self):
        try:
            return self.run.runstatbarc.first().sample_time
            # return self.RunDetails.last().minKNOW_start_time
        except AttributeError:
            return "undefined"

    def last_read(self):
        try:
            return self.run.runstatbarc.last().sample_time
            # return self.reads.last().created_date
        except AttributeError:
            # return "undefined"
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


class RunStatisticBarcode(models.Model):

    run = models.ForeignKey(
        Run,
        on_delete=models.CASCADE,
        related_name='runstatbarc',
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

    is_pass = models.BooleanField(
        default=True
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

    quality_sum = models.DecimalField(
        decimal_places=2,
        max_digits=12,
        default=0
    )

    channel_presence = models.CharField(
        max_length=3000,
        default='0' * 3000
    )

    class Meta:
        verbose_name = 'Run Statistics Barcode'
        verbose_name_plural = 'Run Statistics Barcodes'
        db_table = 'run_statistics_barcode'

    def __str__(self):
        return "{} {} {} {}".format(
            self.run,
            self.sample_time,
            self.type,
            self.barcode
        )

    def number_active_channels(self):
        return len(self.channel_presence.replace('0', ''))


class ChannelSummary(models.Model):

    run = models.ForeignKey(
        Run,
        on_delete=models.CASCADE,
        related_name='runchansum'
    )

    channel = models.IntegerField(

    )

    read_count = models.BigIntegerField(
        default=0
    )

    read_length = models.BigIntegerField(
        default=0
    )

    class Meta:
        db_table = 'channel_summary'

    def __str__(self):
        return "{} {} {}".format(self.run, self.channel, self.read_count)


class HistogramSummary(models.Model):

    BIN_WIDTH = 900

    run = models.ForeignKey(
        Run,
        on_delete=models.CASCADE
    )

    barcode = models.ForeignKey(
        Barcode,
        on_delete=models.CASCADE
    )

    read_type = models.ForeignKey(
        FastqReadType
    )

    is_pass = models.BooleanField(
        default=True
    )

    bin_index = models.BigIntegerField(

    )

    read_count = models.BigIntegerField(
        default=0
    )

    read_length = models.BigIntegerField(
        default=0
    )

    class Meta:
        db_table = 'histogram_summary'

    def __str__(self):
        return "{} {} {}".format(self.run, self.read_type, self.bin_index)


class RunSummaryBarcode(models.Model):

    run = models.ForeignKey(
        Run,
        on_delete=models.CASCADE,
        related_name='runsummariesbarcodes'
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

    is_pass = models.BooleanField(
        default=True
    )  # pass = true, fail = false

    quality_sum = models.DecimalField(
        decimal_places=2,
        max_digits=12,
        default=0
    )

    read_count = models.IntegerField(
        default=0
    )

    total_length = models.BigIntegerField(
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
        db_table = 'run_summary_barcode'

    def __str__(self):
        return "{} {} {} {} {}".format(
            self.run,
            self.total_length,
            self.read_count,
            self.type,
            self.barcode
        )

    def number_active_channels(self):
        return len(self.channel_presence.replace('0', ''))


@receiver(post_save, sender=Run)
def create_run_barcodes(sender, instance=None, created=False, **kwargs):

    if created:

        Barcode.objects.update_or_create(
            run=instance,
            name='All reads'
        )

        Barcode.objects.update_or_create(
            run=instance,
            name='No barcode'
        )


@receiver(post_save, sender=GroupRun)
def create_grouprun_barcodes(sender, instance=None, created=False, **kwargs):

    if created:

        GroupBarcode.objects.get_or_create(
            grouprun=instance,
            name='All reads'
        )

        GroupBarcode.objects.get_or_create(
            grouprun=instance,
            name='No barcode'
        )


@receiver(post_save, sender=Barcode)
def add_barcode_to_groupbarcode(sender, instance=None, created=False, **kwargs):

    if created:

        check_barcode_groupbarcode(instance.run)


def run_groupruns_changed(sender, action, instance, **kargs):

    if action == 'post_add':

        print('run_groupruns_changed')

        check_barcode_groupbarcode(instance)


m2m_changed.connect(run_groupruns_changed, sender=Run.groupruns.through)


def check_barcode_groupbarcode(instance):
    """
    This function checks if all barcodes from this run (instance)
    are associated with the groupruns (through groupbarcodes).

    For each barcode of a run, it goes through all groupbarcodes of all groupruns associated with the run
    and checks if there a groupbarcode's name equal to the barcode name.

    (a) If yes, no action is required.
    (b) If not, but the groupbarcode exists (groupbarcode name is equal to barcode name), it add a manytomany entry.
    (c) If not, and the groupbarcode does not exist, it creates a new groupbarcode and add to the grouprun.
    """

    barcode_list = instance.barcodes.all()

    for grouprun in instance.groupruns.all():

        # print('> grouprun: {}'.format(grouprun.name))

        for barcode in barcode_list:

            # print('>> barcode: {}'.format(barcode.name))

            has_association = False

            groupbarcode_list = grouprun.groupbarcodes.all()

            for groupbarcode in groupbarcode_list:

                # print('>>> groupbarcode: {}'.format(groupbarcode.name))

                if barcode.name == groupbarcode.name:

                    # print('equal names')

                    if barcode not in groupbarcode.barcodes.all():

                        barcode.groupbarcodes.add(groupbarcode)  # case (b)

                        # print('add barcode to groupbarcode')

                    has_association = True

            if not has_association:

                new_groupbarcode = GroupBarcode.objects.create(
                    grouprun=grouprun,
                    name=barcode.name
                )

                barcode.groupbarcodes.add(new_groupbarcode)  # case (c)

                # print('created new groupbarcode and added barcode')
