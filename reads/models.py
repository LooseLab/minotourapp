from django.conf import settings
from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver
from rest_framework.authtoken.models import Token

from minknow_data.models import Minion, Run, Flowcell
from reference.models import ReferenceInfo


class UserOptions(models.Model):
    """
    :purpose: Store information about the user options

    Fields:

    :owner: (django user model object)
    :twitterhandle: (str)
    :tweet: (bool)
    :email: (bool)
    """

    owner = models.OneToOneField(

        settings.AUTH_USER_MODEL,
        related_name='extendedopts',
        on_delete=models.CASCADE,
    )

    twitterhandle = models.CharField(

        max_length=64
    )

    tweet = models.BooleanField(

        default=False
    )

    email = models.BooleanField(

        default=False
    )

    def __str__(self):
        return "{}".format(str(self.owner))


class FastqFile(models.Model):
    """
    :purpose: Provide a lookup to see if a fastq file has been seen before.

    Fields:

    :fastqname: Name of the fastq file.
    :runid: unique run id generated by MinKNOW.
    :checksum: the file checksum as recorded by MinFQ
    """
    name = models.CharField(

        max_length=256
    )

    run = models.ForeignKey(
        Run,
        blank=True,
        null=True,
        on_delete = models.CASCADE,
        related_name = 'fastqfiles'
    )

    runid = models.CharField(
        max_length=64
    )

    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        related_name='fastqfiles',
        on_delete=models.DO_NOTHING,

    )

    md5 = models.CharField(
        blank=True,
        null=True,
        max_length=256
    )

    class Meta:
        unique_together = (("name", "run"),)

    def __str__(self):
        return "{} {}".format(str(self.name),str(self.run))


class Barcode(models.Model):
    """
    :purpose: To store barcode names associated with a particular run.

    Fields:

    :run: Foreign key linking a barcode to a specific run.
    :name: Barcode Name - derived from the read name.

    """

    run = models.ForeignKey(

        Run,
        on_delete=models.CASCADE,
        related_name='barcodes'
    )

    name = models.CharField(

        max_length=32
    )

    def __str__(self):
        return "{} {} {}".format(self.run, self.run.runid, self.name)


class FastqReadType(models.Model):
    """
    :purpose: Store information about the FASTQ read type (eg template/complemnt/1d^2)

    Fields:

    :name: (str) FASTQ read type
    """

    name = models.CharField(

        max_length=16
    )

    class Meta:
        verbose_name = 'FASTQ Read Type'
        verbose_name_plural = 'FASTQ Read Types'

    def __str__(self):
        return self.name


class FastqRead(models.Model):
    """
    :purpose: Each read has a fastqread object. Contains the header information broken down, and some metadata about the read.

    Fields:
    :run: Foreign key linked to run object
    :read_id: The read ID
    :read: Read number auto incrementing in fastq header
    :channel: Channel number that the read was sequenced on
    :barcode: The barcode identifer/number if any
    :sequence_length:
    :quality_average:
    :is_pass: Whether the read passed QC or not
    :type: FK to fastqreadtype, example 1d^2
    :start_time:
    :created_date:
    :modified_Date:
    """
    run = models.ForeignKey(
        Run,
        on_delete=models.CASCADE,
        related_name='reads'
    )

    read_id = models.CharField(
        db_index=True,
        max_length=96
    )

    read = models.IntegerField(

        db_index=True
    )

    channel = models.IntegerField(

        db_index=True
    )

    barcode = models.ForeignKey(
        Barcode,
        # on_delete=models.CASCADE,
        related_name='reads',
        null=True,
        on_delete=models.DO_NOTHING,
    )

    rejected_barcode = models.ForeignKey(
        Barcode,
        on_delete=models.DO_NOTHING,
        related_name="rejection_status",
        null=True,
        blank=True,
    )

    barcode_name = models.CharField(
        db_index=True,
        max_length=32,
    )

    sequence_length = models.BigIntegerField(
        null=True,
        blank=True,
        db_index=True
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
        FastqReadType,
        on_delete=models.DO_NOTHING,
    )

    start_time = models.DateTimeField(
        db_index=True
    )

    created_date = models.DateTimeField(
        auto_now_add=True
    )

    modified_date = models.DateTimeField(
        auto_now=True
    )


    fastqfile = models.ForeignKey(
        FastqFile,
        on_delete=models.CASCADE,
        related_name='reads',
        null=True,
    )

    sequence = models.TextField(
        blank=True,
        null=True
    )
    
    quality = models.TextField(
        blank=True,
        null=True,
    )

    flowcell = models.ForeignKey(
        Flowcell,
        on_delete=models.CASCADE,
        null=True,
        blank=True,
    )

    class Meta:
        verbose_name = 'FASTQ Read'
        verbose_name_plural = 'FASTQ Read'

    def __str__(self):
        #return "1"
        return str(self.read_id)


class RunStatisticBarcode(models.Model):  # TODO to be removed

    run = models.ForeignKey(
        Run,
        on_delete=models.CASCADE,
        related_name='runstatbarc',
        null=True
    )

    type = models.ForeignKey(
        FastqReadType,
        on_delete=models.DO_NOTHING,
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

    channel_count = models.IntegerField(
        default=0
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


class FlowcellStatisticBarcode(models.Model):
    """
    model to contain the flowcell statistic barcodes data for the chancalc page
    """
    flowcell = models.ForeignKey(

        Flowcell,
        on_delete=models.CASCADE,
        related_name='statistics',
        null=True
    )

    read_type_name = models.CharField(

        max_length=32
    )

    barcode_name = models.CharField(

        max_length=32
    )

    rejection_status = models.CharField(
        max_length=32,
        default="Sequenced",
        null=True,
        blank=True,
    )

    status = models.CharField(

        max_length=32
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

    channel_count = models.IntegerField(
        default=0
    )

    class Meta:
        verbose_name = 'Flowcell Statistics Barcode'
        verbose_name_plural = 'Flowcell Statistics Barcodes'
        db_table = 'flowcell_statistics_barcode'
        constraints = [
            models.UniqueConstraint(
                fields=['flowcell', 'sample_time','barcode_name', 'rejection_status', 'read_type_name', 'status'],
                name='flowcell_stats_barc_summ')
        ]
        #unique_together = ('flowcell','sample_time','barcode_name','rejection_status','read_type_name','status')

    def __str__(self):
        return "{} {} {} {}".format(
            self.flowcell.name,
            self.sample_time,
            self.read_type_name,
            self.barcode_name
        )

    def number_active_channels(self):
        return len(self.channel_presence.replace('0', ''))


class ChannelSummary(models.Model):  # don't document

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


class FlowcellChannelSummary(models.Model):  # TODO to be deleted
    """
    The flowcell channel model for the channel visualistions on te chancalc page
    """
    flowcell = models.ForeignKey(
        Flowcell,
        on_delete=models.CASCADE,
        related_name='channel_summaries'
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
        db_table = 'flowcell_channel_summary'
        constraints = [
            models.UniqueConstraint(fields=['flowcell', 'channel'], name='flowcell_channel')
        ]

    def __str__(self):
        return "{} {} {}".format(self.flowcell, self.channel, self.read_count)


class HistogramSummary(models.Model):

    BIN_WIDTH = 1000

    run = models.ForeignKey(
        Run,
        on_delete=models.CASCADE
    )

    barcode = models.ForeignKey(
        Barcode,
        on_delete=models.CASCADE
    )

    read_type = models.ForeignKey(
        FastqReadType,
        on_delete=models.DO_NOTHING,
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


class FlowcellHistogramSummary(models.Model):
    """
    :purpose: Summarise number of reads (read_count) and read length (read_length) in bins of 900 bases width.

    Fields:

    :flowcell: (Flowcell) Foreign key to Flowcell
    :read_type_name: (FastReadType) FastqReadType name
    :barcode_name: (Barcode) Barcode name
    :status: (boolean) Pass or Fail; originates from FastqRead is_pass attribute
    :bin_index: (int) Bin position
    :read_count: (float) Total number of all reads from all runs of the flowcell
    :read_length: (float) Sum of the length of all reads from all runs of the flowcell
    """

    BIN_WIDTH = 1000

    flowcell = models.ForeignKey(

        Flowcell,
        on_delete=models.CASCADE
    )

    read_type_name = models.CharField(

        max_length=32
    )

    barcode_name = models.CharField(

        max_length=32
    )

    rejection_status = models.CharField(
        max_length=32,
        default="Sequenced",
        null=True,
        blank=True,
    )

    status = models.CharField(

        max_length=32
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
        db_table = 'flowcell_histogram_summary'
        constraints = [
            models.UniqueConstraint(fields=['flowcell','barcode_name','rejection_status','read_type_name','status','bin_index'], name='flowcell_hist_summ')
        ]
        #unique_together =('flowcell','barcode_name','rejection_status','read_type_name','status','bin_index')

    def __str__(self):
        return "{} {} {}".format(
            self.flowcell,
            self.read_type_name,
            self.status,
            self.bin_index
        )


class RunSummary(models.Model):
    """
    This class keep summary information about individual runs
    There should be a task updating this information regularly
    """

    run = models.OneToOneField(

        Run,
        on_delete=models.CASCADE,
        related_name='summary'
    )

    read_count = models.IntegerField(

        null=True,
        blank=True
    )

    total_read_length = models.BigIntegerField(

        null=True,
        blank=True
    )

    max_read_length = models.BigIntegerField(

        null=True,
        blank=True
    )

    avg_read_length = models.BigIntegerField(

        null=True,
        blank=True
    )

    min_read_length = models.BigIntegerField(

        null=True,
        blank=True
    )

    first_read_start_time = models.DateTimeField(

        null=True,
        blank=True
    )

    last_read_start_time = models.DateTimeField(

        null=True,
        blank=True
    )

    last_read = models.BigIntegerField(

        default=0
    )

    running = models.BooleanField(
        default=False
    )



    class Meta:

        verbose_name = 'Run Summary'
        verbose_name_plural = 'Run Summaries'
        db_table = 'run_summary'

    def __str__(self):

        return self.run.name


class RunSummaryBarcode(models.Model):  # TODO to be deleted

    run = models.ForeignKey(
        Run,
        on_delete=models.CASCADE,
        related_name='runsummariesbarcodes'
    )

    type = models.ForeignKey(
        FastqReadType,
        on_delete=models.DO_NOTHING,
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

    channel_count = models.IntegerField(
        default=0
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
        """
        TODO Remove
        """

        return len(self.channel_presence.replace('0', ''))


class FlowcellSummaryBarcode(models.Model):
    """
    :purpose: Summarise information from runs by flowcell, barcode name, fastq read type, and status (pass or fail).
     There is one record per flowcell. Most of the charts in the flowcell page use this data.

    Fields:

    :flowcell: (Flowcell) Foreign key to Flowcell
    :read_type_name: (FastReadType) FastqReadType name
    :barcode_name: (Barcode) Barcode name
    ;rejection_status: (Barcode) Foreign key to whether the read was accepted or unblocked
    :status: (boolean) Pass or Fail; originates from FastqRead is_pass attribute
    :quality_sum: (float) # TODO
    :read_count: (float) Total number of all reads from all runs of the flowcell
    :total_length: (float) Sum of the length of all reads from all runs of the flowcell
    :max_length: (float) Maximum read length of the flowcell
    :min_length: (float) Minimum read length of the flowcell
    :channel_presence: (str) Sequence of 3000 zeros and ones representing the presence or absence of a strand
    :channel_count: (int) Retuns the sum of ones in the channel_presence
    """

    flowcell = models.ForeignKey(

        Flowcell,
        on_delete=models.CASCADE,
        related_name='flowcellsummariesbarcodes'
    )

    read_type_name = models.CharField(

        max_length=32
    )

    barcode_name = models.CharField(

        max_length=32
    )

    rejection_status = models.CharField(
        max_length=32,
        default="Sequenced",
        null=True,
        blank=True,
    )

    status = models.CharField(

        max_length=32
    )

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

    channel_count = models.IntegerField(

        default=0
    )

    class Meta:
        verbose_name = 'Flowcell Summary'
        verbose_name_plural = 'Flowcell Summary'
        db_table = 'flowcell_summary_barcode'
        constraints = [
            models.UniqueConstraint(
                fields=['flowcell', 'barcode_name', 'rejection_status', 'read_type_name', 'status'],
                name='flowcell_summ')
        ]
        #unique_together = ('flowcell','barcode_name','rejection_status','read_type_name','status')

    def __str__(self):
        return "{} {} {} {} {}".format(
            self.flowcell,
            self.total_length,
            self.read_count,
            self.read_type_name,
            self.barcode_name
        )

    def average_read_length(self):
        """
        TODO
        """

        if self.read_count > 0:
            return self.total_length / self.read_count
        else:
            return 0

    def number_active_channels(self):
        """
        TODO
        """

        return len(self.channel_presence.replace('0', ''))


class JobType(models.Model):

    name = models.CharField(

        max_length=256,
    )

    description = models.TextField(

        max_length=256,
        blank=True,
        null=True,
    )

    long_description = models.TextField(

        blank=True,
        null=True,
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

        default=True,
    )

    def __str__(self):

        return "{}".format(self.name)


class JobMaster(models.Model):
    """
    The JobMaster model is used to store the tasks that we wish celery to pick up and run
    """
    run = models.ForeignKey(
        Run,
        on_delete=models.DO_NOTHING,
        related_name='run_jobs',
        null=True,
        blank=True,
    )
    flowcell = models.ForeignKey(
        Flowcell,
        on_delete=models.CASCADE,
        related_name='flowcell_jobs',
        null=True,
        blank=True,
    )
    job_type = models.ForeignKey(
        JobType,
        on_delete=models.DO_NOTHING,
        related_name='task_name',
    )
    reference = models.ForeignKey(
        ReferenceInfo,
        on_delete=models.SET_NULL,
        related_name='reference_job',
        null=True,
        blank=True,
    )
    barcode = models.ForeignKey(
        Barcode,
        on_delete=models.CASCADE,
        related_name="jobs_barcode",
        null=True,
        blank=True
    )
    last_read = models.BigIntegerField(
        default=0,
    )
    tempfile_name = models.CharField(
        max_length=256,
        blank=True,
        null=True,
    )
    read_count = models.BigIntegerField(
        default=0
    )
    complete = models.BooleanField(
        default=False,
    )
    running = models.BooleanField(
        default=False,
    )
    paused = models.BooleanField(
        default=False
    )
    from_database = models.BooleanField(
        default=False
    )
    iteration_count = models.IntegerField(
        null=True,
        default=0,
    )
    target_set = models.CharField(
        default=None,
        null=True,
        max_length=100,
        blank=True,
    )
    def __str__(self):
        if self.run is not None:
            return f"Run: {self.run} JobType: {self.job_type} Id: {self.id} Barcode: {self.barcode}"
        return f"Flowcell: {self.flowcell_id} JobType: {self.job_type} Id: {self.id} Barcode: {self.barcode}"


class SampleTag(models.Model):

    flowcell = models.ForeignKey(

        Flowcell,
        on_delete=models.CASCADE,
        related_name='samples',
    )

    key_text = models.CharField(
        max_length=64
    )

    value_text = models.TextField(

    )


@receiver(post_save, sender=settings.AUTH_USER_MODEL)
def create_auth_token(sender, instance=None, created=False, **kwargs):
    """
    TODO
    """

    if created:
        Token.objects.create(user=instance)


@receiver(post_save, sender=Run)
def create_run_barcodes(sender, instance=None, created=False, **kwargs):
    """
    TODO No barcode shouldn't be part of default right?? Come back to this, I think I'm right
    """

    if created:
        Barcode.objects.update_or_create(
            run=instance,
            name='All reads'
        )

        Barcode.objects.update_or_create(
            run=instance,
            name='No barcode'
        )

