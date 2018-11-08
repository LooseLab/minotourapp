import datetime

import pytz
from django.conf import settings
from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver
from rest_framework.authtoken.models import Token


class Experiment(models.Model):
    """
    Define a group of runs
    """

    name = models.CharField(

        max_length=32
    )

    owner = models.ForeignKey(

        settings.AUTH_USER_MODEL,
        related_name='experiments'
    )

    created_date = models.DateTimeField(

        auto_now_add=True
    )

    modified_date = models.DateTimeField(

        auto_now=True
    )


class Flowcell(models.Model):

    name = models.CharField(
        max_length=256,
    )

    sample_name = models.CharField(
        max_length=256,
    )

    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        related_name='flowcells'
    )

    size = models.IntegerField(
        default=512
    )

    is_active = models.BooleanField(
        default=True
    )

    start_time = models.DateTimeField(
        auto_now_add=True
    )

    number_reads = models.BigIntegerField(
        default=0
    )

    number_reads_processed = models.BigIntegerField(
        default=0
    )

    number_runs = models.IntegerField(
        default=0
    )

    number_barcodes = models.IntegerField(
        default=0
    )

    average_read_length = models.BigIntegerField(
        default=0
    )

    total_read_length = models.BigIntegerField(
        default=0
    )

    created_date = models.DateTimeField(
        auto_now_add=True
    )

    modified_date = models.DateTimeField(
        auto_now=True
    )

    experiments = models.ManyToManyField(
        Experiment,
        blank=True
    )

    has_fastq = models.BooleanField(

        default=True
    )

    def barcodes(self):

        barcode_set = set()

        for run in self.runs.all():

            for barcode in run.barcodes.all():
                barcode_set.add(barcode)

        return barcode_set

    def __str__(self):
        return "{} {}".format(self.name, self.id)


class MinION(models.Model):
    minION_name = models.CharField(
        max_length=64
    )

    name = models.CharField(
        max_length=64,
        blank=True,
        null=True
    )

    owner = models.ForeignKey(

        settings.AUTH_USER_MODEL,
        related_name='minIONs'
    )

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
            # return reverse('minIONrunstats_list', args=[self.minionrun.last().id])
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


class GroupRun(models.Model):  # TODO don't document

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

        max_length=32
    )

    device = models.CharField(

        max_length=10,
        choices=NANOPORE_DEVICES,
        default=MINION
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
    """
    :purpose: Store a new job for the minion to pickup and perform  TODO could probably be expanded with example job, why it exists

    Fields:
    :minion: The minion name as stored in the Minion table
    :owner: The Owner as the user
    :job: The Job that you want the minion to perform
    :custom: TODO to be completed
    :complete: Whether the job is complete or not
    :created_date: The date that the job was created
    :modified_date: The date that this as last modified
    """
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
        related_name='extendedopts'
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


class Run(models.Model):
    """
    :purpose: Store information about each Run generated during the sequencing process in a flowcell. A flowcell can contain multiple runs, usually alternating short ones (muxscan) with long ones (real sequence data). 

    Fields:

    :minion: the MinION that contains the flowcell used the produce the run data
    :flowcell: the flowcell used to produce the run data
    :owner: the run's owner; usually provides the token used in the minotour client
    :groupruns: to be removed
    :name: defines the run's name; this information usually come from the -k parameter in the minotour client
    :runid: unique identified generated by MinKnow
    :is_barcoded: true if the run has barcoded read, false otherwise
    :has_fastq: true if the sequence and quality were sent by the minotour client (stored at FastqReadExtra)
    :active: true if the Run is active; defined by a set of criteria and updated by a task
    :to_delete: to be removed
    :start_time: returns the time of the earliest read; updated by web.tasks.update_run_start_time
    """

    minion = models.ForeignKey(

        MinION,
        blank=True,
        null=True,
        related_name='runs'
    )

    flowcell = models.ForeignKey(

        Flowcell,
        related_name='runs',
    )

    owner = models.ForeignKey(

        settings.AUTH_USER_MODEL,
        related_name='runs'
    )

    groupruns = models.ManyToManyField(

        GroupRun,
        blank=True,
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

    start_time = models.DateTimeField(

        blank=True,
        null=True
    )

    class Meta:

        verbose_name = 'Run'
        verbose_name_plural = 'Runs'

    def __str__(self):

        return "{} - {}".format(self.name, self.runid)

    def flowcell_name(self):
        """
        Returns flowcell name - to be removed and replaces by run's instance.flowcell.name
        """

        try:
            return self.flowcell.name

        except AttributeError:
            return "undefined"

    def last_entry(self):
        """
        Returns the time of the most recent RunStat - to be removed and replaces by a field and updated by a task
        """

        try:
            return self.RunStats.last().created_date
        except AttributeError:
            # return "undefined"
            olddate = datetime.datetime(1, 1, 1, 1, 1, 1, 1, pytz.utc)
            return olddate

    def last_read(self):
        """
        Returns the time of the most recent read - to be removed and replaces by a field and updated by a task
        """

        try:
            return self.reads.last().created_date
        except AttributeError:
            # return "undefined"
            olddate = datetime.datetime(1, 1, 1, 1, 1, 1, 1, pytz.utc)
            return olddate

    def sample_name(self):
        """
        TODO Matt - what is the connection between Run and RunDetails.
        """

        try:
            return self.RunDetails.last().minKNOW_sample_name
        except AttributeError:
            return "undefined"

    def minKNOW_flow_cell_id(self):
        """
        TODO Matt - what is the connection between Run and RunDetails. Runs belong to flowcells at the moment.
        """
        try:
            return self.RunDetails.last().minKNOW_flow_cell_id
        except AttributeError:
            return "undefined"

    def minKNOW_version(self):
        """
        TODO Matt
        """
        try:
            return self.RunDetails.last().minKNOW_version
        except AttributeError:
            return "undefined"

    def max_channel(self):
        """
        TODO Matt - Maybe removed and replace by field and task
        """
        try:
            return self.runchansum.order_by('channel').last().channel
        except AttributeError:
            return "undefined"

    def flowcell_type(self):
        """
        TODO Matt - Maybe removed and replace by field and task
        """
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

    def computer_name(self):
        """
        TODO Matt - Maybe removed and replace by field and task
        """

        if self.minion and self.minion.currentrundetails and self.minion.currentrundetails.last():

            return self.minion.currentrundetails.last().minKNOW_computer

        else:

            return None


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
        related_name='fastqfiles'
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
    :barcodegroup: Foreign key linking to a specific barcode group. #ToDo Remove this field.
    :groupbarcodes: Foreign key linking to a group barcode. #ToDo Remove this field.
    :name: Barcode Name - derived from the read name.

    :run: Returns the run that owns the barcode.
    :barcodegroup: To be removed.
    :groupbarcodes: To be removed.
    :name: returns the name as identified in the header of the fastq file.
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


class MinIONStatus(models.Model):  # TODO Rename class for logical consistency
    """
    :purpose: To store data on the current status of a specific MinION at a given moment in time.
              This includes detailed pon the connected computer, the space available and other key run details.

    Fields:

    :minION: Foreign key linking to a MinION
    :minKNOW_status: Current status of the MinION device (Ready, Starting, Processing or Finishing) as reported by MinKNOW
    :minKNOW_current_script: Current or last run sequencing script from MinKNOW for this MinION
    :minKNOW_sample_name: Current or last used Sample Name
    :minKNOW_exp_script_purpose: Experiment purpose as identified by MinKNOW
    :minKNOW_flow_cell_id: The unique identifier for a flow cell
    :minKNOW_run_name: The run name as determined by MinKNOW
    :minKNOW_hash_run_id: A hash of the same run name
    :minKNOW_script_run_id: A script identifier. The format of which is output by MinKNOW and variable.
    :minKNOW_real_sample_rate: The rate at which data have been sampled for this run
    :minKNOW_asic_id: The specific asic identifier
    :minKNOW_total_drive_space: The total hard drive space available to the currently sequencing computer.
    :minKNOW_disk_space_till_shutdown: How much space is left till the sequencer shuts down.
    :minKNOW_disk_available: The total amount of space left on the device.
    :minKNOW_warnings: If minKNOW is about to shutdown a warning is added here. #todo Contact the user if the computer is about to shut down.

    """

    minION = models.OneToOneField(

        MinION,
        related_name='currentdetails'

    )

    minKNOW_status = models.CharField(

        max_length=64

    )

    minKNOW_current_script = models.CharField(

        max_length=256,
        blank=True,
        null=True

    )

    minKNOW_sample_name = models.CharField(

        max_length=256,
        blank=True,
        null=True

    )

    minKNOW_exp_script_purpose = models.CharField(

        max_length=256,
        blank=True,
        null=True

    )

    minKNOW_flow_cell_id = models.CharField(

        max_length=64,
        blank=True,
        null=True

    )

    minKNOW_run_name = models.CharField(

        max_length=256,
        blank=True,
        null=True

    )

    minKNOW_hash_run_id = models.CharField(

        max_length=256,
        blank=True,
        null=True

    )

    minKNOW_script_run_id = models.CharField(

        max_length=256,
        blank=True,
        null=True

    )

    minKNOW_real_sample_rate = models.IntegerField(
        blank=True,
        null=True

    )

    # minKNOW_voltage_offset = models.IntegerField(blank=True, null=True)

    # minKNOW_yield = models.IntegerField(blank=True, null=True)

    minKNOW_asic_id = models.CharField(

        max_length=256,
        blank=True,
        null=True

    )

    minKNOW_total_drive_space = models.FloatField(

        blank=True,
        null=True

    )

    minKNOW_disk_space_till_shutdown = models.FloatField(

        blank=True,
        null=True

    )

    minKNOW_disk_available = models.FloatField(

        blank=True,
        null=True

    )

    minKNOW_warnings = models.BooleanField(

        default=False

    )

    class Meta:
        verbose_name = 'MinION Status'
        verbose_name_plural = 'MinION Status'

    def __str__(self):
        return "{} {}".format(self.minION, self.minKNOW_status)


class MinIONRunStats(models.Model):  # Todo consider merging in to one object with
    """
    :purpose: ???

    Fields:

    :minion: Foreign key to link to a specific MinION run.
    :run_id: Foreign key to link to a specific run id.
    :sample_time: Date Time field identifying the moment of data collection.
    :event_yield: Sequencing yield as approximated in events by MinKNOW.
    :asic_temp: Current temperature of the ASIC.
    :heat_sink_temp: Current temperature of the heat sink.
    :voltage_value: Current set voltage on the device.
    :mean_ratio: A ratio of currents for in strand versus open pore. #ToDO not currently reported.
    :open_pore: Count of pore classification type (as named)
    :in_strand: Count of pore classification type (as named)
    :multiple: Count of pore classification type (as named)
    :unavailable: Count of pore classification type (as named)
    :unknown: Count of pore classification type (as named)
    :adapater: Count of pore classification type (as named)
    :pending_mux_change: Count of pore classification type (as named)
    :unclassified: Count of pore classification type (as named)
    :below: Count of pore classification type (as named)
    :unblocking: Count of pore classification type (as named)
    :above: Count of pore classification type (as named)
    :good_single: Count of pore classification type (as named)
    :saturated: Count of pore classification type (as named)
    :inrange: Count of pore classification type (as named)
    :strand: Count of pore classification type (as named)
    :minKNOW_read_count: Count of reads seen at this point.
    :minKNOW_histogram_values: Reporting the values of the histogram from MinKNOW
    :minKNOW_histogram_bin_width: Measure of the bin width for the above histogram.
    :created_date: Created Date.
    """

    minION = models.ForeignKey(

        MinION,
        related_name='currentrunstats'
    )

    run_id = models.ForeignKey(

        Run,
        related_name='RunStats'
    )

    sample_time = models.DateTimeField(


    )

    event_yield = models.BigIntegerField(

        default=0
    )

    asic_temp = models.FloatField(

        default=0
    )

    heat_sink_temp = models.FloatField(

        default=0
    )

    voltage_value = models.FloatField(

        default=0
    )

    mean_ratio = models.FloatField(

        default=0
    )

    open_pore = models.FloatField(

        default=0
    )

    in_strand = models.FloatField(

        default=0
    )

    multiple = models.IntegerField(

        default=0
    )

    unavailable = models.IntegerField(

        default=0
    )

    unknown = models.IntegerField(

        default=0
    )

    adapter = models.IntegerField(

        default=0
    )

    pending_mux_change = models.IntegerField(

        default=0
    )

    unclassified = models.IntegerField(

        default=0
    )

    below = models.IntegerField(

        default=0
    )

    unblocking = models.IntegerField(

        default=0
    )

    above = models.IntegerField(

        default=0
    )

    good_single = models.IntegerField(

        default=0
    )

    saturated = models.IntegerField(

        default=0
    )

    inrange = models.IntegerField(

        default=0
    )

    strand = models.IntegerField(

        default=0
    )

    minKNOW_read_count = models.IntegerField(

        default=0
    )

    minKNOW_histogram_values = models.TextField(

        blank=True, null=True
    )

    minKNOW_histogram_bin_width = models.IntegerField(

        default=900
    )

    created_date = models.DateTimeField(

        auto_now_add=True
    )

    class Meta:
        verbose_name = 'MinION Run Stats'
        verbose_name_plural = 'MinION Run Stats'

    def __str__(self):
        return "{} {} {}".format(self.minION, self.run_id, self.sample_time)

    ## This is something to look at for optimisation
    def occupancy(self):
        if (self.strand > 0 and self.inrange > 0):
            occupancy = round(((self.strand + self.adapter) / (self.strand + self.adapter + self.good_single)) * 100)
        else:
            occupancy = 0
        return occupancy


class MinIONRunStatus(models.Model):
    """
    :purpose: To store data on the current status of a specific MinION at a given moment in time.
              This includes detailed pon the connected computer, the space available and other key run details.

    Fields:

    :minION:
    :minKNOW_current_script:
    :minKNOW_sample_name:
    :minKNOW_exp_script_purpose:
    :minKNOW_flow_cell_id:
    :minKNOW_version:
    :minKNOW_run_name:
    :run_id:
    :minKNOW_hash_run_id:
    :minKNOW_script_run_id:
    :minKNOW_real_sample_rate:
    :minKNOW_asic_id:
    :minKNOW_start_time:
    :minKNOW_colours_string:
    :minKNOW_computer:
    :experiment_type:
    :experiment_id:
    :fast5_output_fastq_in_hdf:
    :fast5_raw:
    :fast5_reads_per_folder:
    :fastq_enabled:
    :fastq_reads_per_file:
    :filename:
    :flowcell_type:
    :kit_classification:
    :local_basecalling:
    :sample_frequency:
    :sequencing_kit:
    :user_filename_input:
    """

    minION = models.ForeignKey(

        MinION,
        related_name='currentrundetails'
    )

    minKNOW_current_script = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    minKNOW_sample_name = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    minKNOW_exp_script_purpose = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    minKNOW_flow_cell_id = models.CharField(

        max_length=64,
        blank=True,
        null=True
    )

    minKNOW_version = models.CharField(

        max_length=64,
        blank=True,
        null=True
    )

    minKNOW_run_name = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    run_id = models.ForeignKey(

        Run,
        related_name='RunDetails'
    )

    minKNOW_hash_run_id = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    minKNOW_script_run_id = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    minKNOW_real_sample_rate = models.IntegerField(

        blank=True,
        null=True
    )

    minKNOW_asic_id = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    minKNOW_start_time = models.DateTimeField(

        blank=True,
        null=True
    )

    minKNOW_colours_string = models.TextField(

        blank=True,
        null=True
    )

    minKNOW_computer = models.TextField(

        max_length=128,
        blank=True,
        null=True
    )

    experiment_type = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    experiment_id = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    fast5_output_fastq_in_hdf = models.IntegerField(

        blank=True,
        null=True
    )

    fast5_raw = models.IntegerField(

        blank=True,
        null=True
    )

    fast5_reads_per_folder = models.IntegerField(

        blank=True,
        null=True
    )

    fastq_enabled = models.IntegerField(

        blank=True,
        null=True
    )

    fastq_reads_per_file = models.IntegerField(

        blank=True,
        null=True
    )

    filename = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    flowcell_type = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    kit_classification = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    local_basecalling = models.IntegerField(

        blank=True,
        null=True
    )

    sample_frequency = models.IntegerField(

        blank=True,
        null=True
    )

    sequencing_kit = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    user_filename_input = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    class Meta:
        unique_together = ("minION", "minKNOW_hash_run_id", "minKNOW_exp_script_purpose")
        verbose_name = 'MinION Run Status'
        verbose_name_plural = 'MinION Run Status'

    def __str__(self):
        return "{} {} {}".format(self.minION, self.minKNOW_current_script, self.run_id)

    def minION_name(self):
        return self.minION.minION_name


class MinIONEventType(models.Model):
    """
    :purpose: Provides a list of states that a MinION can be in at any moment in time.

    :Fields:

    :name: The name of the event.
    """

    name = models.CharField(

        max_length=64
    )

    class Meta:
        verbose_name = 'MinION Event Type'
        verbose_name_plural = 'MinION Event Types'

    def __str__(self):
        return self.name


class MinIONEvent(models.Model):
    """
    :purpose: Report the specific history of events occurring to a minION.

    Fields:

    :computer_name: The name of the computer that the MinION is connected to when an event occurs
    :minION: Foreign key to the MinION id
    :event: Foreign key to the event type
    :datetime: datetime that the event occured
    :created_date: date that this event was first created.
    :modified_date: last time a change was made to this entry.

    """

    computer_name = models.CharField(

        max_length=256
    )

    minION = models.ForeignKey(

        MinION, related_name='events'
    )

    event = models.ForeignKey(

        MinIONEventType
    )

    datetime = models.DateTimeField(

    )

    created_date = models.DateTimeField(

        auto_now_add=True
    )

    modified_date = models.DateTimeField(

        auto_now=True
    )  # These last two fields added to enable auto cleanup of event status for a minION incase of disconnect of client.

    class Meta:
        verbose_name = 'MinION Event'
        verbose_name_plural = 'MinION Events'

    def __str__(self):
        return "{} {} {} {}".format(self.computer_name, self.minION, self.event, self.datetime)


class MinIONScripts(models.Model):
    """
    :purpose: Collect all scripts that are available to run in the used version of minKnow and store them in
    the database

    Fields:
    :minION: A FK linking to the minion entry in the database
    :identifier: TODO Matt what is identifier
    :name: The name of the script
    :experiment_type: TODO  MAtt
    :base_calling: Bool field identifies if script is a base calling or not
    :flow_cell: The flowcell name
    :kit: TODO What is kit
    :experiment_time: TODO MATT
    :event_ratio: TODO Matt
    :kit_category: TODO Matt

    """
    minION = models.ForeignKey(

        MinION,
        related_name='scripts'
    )

    identifier = models.CharField(

        max_length=256
    )

    name = models.CharField(

        max_length=256
    )

    experiment_type = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    base_calling = models.NullBooleanField(

        blank=True,
        null=True
    )

    flow_cell = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )

    kit = models.CharField(

        max_length=256,
        blank=True,
        null=True
    )
    experiment_time = models.IntegerField(

        blank=True,
        null=True
    )
    event_ratio = models.FloatField(

        blank=True,
        null=True
    )

    kit_category = models.CharField(
        max_length=1024,
        blank=True,
        null=True
    )

    class Meta:
        verbose_name = 'MinION Script'
        verbose_name_plural = 'MinION Scripts'

    def __str__(self):
        return "{} {} {}".format(self.minION, self.name, self.identifier)


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

    barcode_name = models.CharField(

        max_length=32
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

    fastqfile = models.ForeignKey(
        FastqFile,
        on_delete=models.CASCADE,
        related_name = 'reads',
        null = True
    )

    class Meta:
        verbose_name = 'FASTQ Read'
        verbose_name_plural = 'FASTQ Read'

    def __str__(self):
        #return "1"
        return str(self.read_id)


class FastqReadExtra(models.Model):
    """
    :purpose: If the user choose to send the read sequence and fastq quality, this is the model that is used to store it

    Fields:
    :fastqread: One to One field linking each fastqread object to its fastqreadextra counterpart
    :sequence: The read Sequence
    :quality: The fastq quality chars
    """
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
        return str(self.id)


class MinionMessage(models.Model):
    """
    :purpose: Store all the Messages that the minion sends over

    Fields:
    :minion: FK linking to the minion details stored in the devices_minion table
    :run: FK linking to the Run detail
    :message: The message things like run started etc.
    :identifier: TODO matt
    :severity: TODO matt
    :timestamp:  The timestamp of when nthemessagewassent
    """
    minion = models.ForeignKey(
        MinION,
        related_name='messages'
    )

    run = models.ForeignKey(
        Run,
        related_name='runmessages',
        blank=True,
        null=True
    )

    message = models.CharField(
        max_length=256,
        blank=True,
        null=True
    )

    identifier = models.CharField(
        max_length=256
    )

    severity = models.CharField(
        max_length=64
    )

    timestamp = models.DateTimeField(

    )

    class Meta:
        # unique_together = ("minion", "run", "timestamp")
        unique_together = ("minion", "message",
                           "timestamp")  # Todo note that we assume the same minION will never generate the same message at the same timestamp
        verbose_name = 'Minion message'
        verbose_name_plural = 'MinIon messages'

    def __str__(self):
        return "{} {} {} {}".format(
            self.minion, self.message, self.severity, self.timestamp)


class RunStatisticBarcode(models.Model):  # TODO to be removed

    run = models.ForeignKey(
        Run,
        on_delete=models.CASCADE,
        related_name='runstatbarc',
        null=True
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

    def __str__(self):
        return "{} {} {}".format(self.flowcell, self.channel, self.read_count)


class HistogramSummary(models.Model):  # don't comment

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

    BIN_WIDTH = 900

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
    :purpose: Summarise information from runs by flowcell, barcode name, fastq read type, and status (pass or fail). There is one record per flowcell. Most of the charts in the flowcell page use this data.

    Fields:

    :flowcell: (Flowcell) Foreign key to Flowcell
    :read_type_name: (FastReadType) FastqReadType name
    :barcode_name: (Barcode) Barcode name
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

        return self.total_length / self.read_count

    def number_active_channels(self):
        """
        TODO
        """

        return len(self.channel_presence.replace('0', ''))


class FlowcellTab(models.Model):

    #
    # tab choices
    #

    SUMMARY = 'SU'
    LIVE_EVENT = 'LE'
    BASECALLED = 'BC'
    READ_DATA = 'RD'
    MAPPING = 'MA'
    ASSEMBLY = 'AS'
    METAGENOMICS = 'ME'
    TASK = 'TA'

    TAB_CHOICES = (
        (SUMMARY, 'Summary'),
        (LIVE_EVENT, 'Live Event Data'),
        (BASECALLED, 'Basecalled Data'),
        (READ_DATA, 'Read Data'),
        (MAPPING, 'Sequence Mapping'),
        (ASSEMBLY, 'Assembly'),
        (METAGENOMICS, 'Metagenomics'),
        (TASK, 'Tasks'),
    )

    flowcell = models.ForeignKey(

        Flowcell,
        related_name='tabs',
        on_delete=models.CASCADE
    )

    tab = models.CharField(

        max_length=2,
        choices=TAB_CHOICES
    )

    def __str__(self):
        return '{} - {}'.format(self.flowcell.id, self.tab)



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
    TODO
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
