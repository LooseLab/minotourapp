import datetime

import pytz
from django.conf import settings
from django.db import models

# Create your models here.
from minotourapp.utils import get_env_variable


class Minion(models.Model):
    """
    Contains info about a Minion Device - such as name, and links together all the information about live data runs
    """
    minION_name = models.CharField(
        max_length=64
    )

    name = models.CharField(
        max_length=64,
        blank=True,
        null=True
    )
    computer = models.CharField(
        max_length=100,
        blank=True,
        null=True
    )
    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        related_name='minIONs',
        on_delete=models.DO_NOTHING,
    )

    class Meta:
        verbose_name = 'MinION'
        verbose_name_plural = 'MinIONs'

    def __str__(self):
        return self.minION_name

    def start_time(self):
        """
        Fetch the run start time according to MinKnow from the linked minionrunstatus model.
        Returns
        -------
        start_time: datetime.datetime
        """
        try:
            start_time = self.currentrundetails.last().minKNOW_start_time
        except AttributeError:
            start_time = "Unknown"

        return start_time

    def status(self):
        try:
            status = self.events.order_by('datetime').last().event.name
        except AttributeError:
            status = "undefined"
        return status

    def last_run(self):
        try:
            return self.currentrundetails.last().id
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

    def minknow_version(self):
        try:
            return self.currentdetails.minknow_version
        except AttributeError as e:
            print(e)
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

    def current_script(self):
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

    def actual_max_val(self):
        """
        Fetch the actual longest read from minKnow by reverse looking up the last MinionRunInfo entry.
        Returns
        -------
        int
        """
        try:
            return self.currentrunstats.order_by('sample_time').last().actual_max_val
        except AttributeError:
            return "Unknown"

    def wells_per_channel(self):
        """
        Fetch the number of wells per channel by looking up the last Minion run info entry.
        Returns
        -------

        """
        try:
            return self.currentrundetails.last().wells_per_channel
        except AttributeError:
            return "Unknown"

    def read_length_type(self):
        """
        Fetch the read length type of the current run on the minion from Minion run info entry.
        Returns
        -------

        """
        try:
            return self.currentrundetails.last().read_length_type
        except AttributeError:
            return "Unknown"

    def target_temp(self):
        """
        Fetch the target temperature of the Minion device.
        Returns
        -------
        int
        """
        try:
            return self.currentrundetails.last().target_temp
        except AttributeError:
            return 0

    def flowcell_type(self):
        """
        Fetch the flowcell product type name
        Returns
        -------
        str

        """
        try:
            return self.currentrundetails.last().flowcell_type
        except AttributeError:
            return "Unknown"

    def experiment_name(self):
        """
        Return the experiment id in string format
        Returns
        -------

        """
        try:
            return self.currentrundetails.last().experiment_id
        except AttributeError:
            return  "Unknown"


class Flowcell(models.Model):
    """
        The model for the django ORM that deals with the creation and retrieval of Flowcells.
        Flowcells are top level organisers under users.
    """
    name = models.CharField(
        max_length=256,
    )

    sample_name = models.CharField(
        max_length=256,
    )

    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        related_name='flowcells',
        on_delete=models.DO_NOTHING,
    )

    max_channel = models.IntegerField(
        default=0
    )

    size = models.IntegerField(
        default=512
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

    has_fastq = models.BooleanField(

        default=True
    )

    last_activity_date = models.DateTimeField(
        auto_now_add=True
    )

    last_activity_message = models.CharField(

        max_length=256,
        default='Created flowcell.'
    )

    archived = models.BooleanField(
        default=False
    )

    sample_description = models.TextField(
        null=True,
        blank=True,
    )

    def barcodes(self):
        """
        Return all the barcodes in runs contained under this flowcell.
        :return:
        """

        barcode_set = set()
        # for each run in all the runs foreign keyed to this flowcell
        for run in self.runs.all():
            # for all the barcodes foreign keyed to this run
            for barcode in run.barcodes.all():
                # add the barcode to the set - guarantees uniqueness
                barcode_set.add(barcode)

        return barcode_set

    def __str__(self):
        """
        Format the string returned when models are inspected in the shell or by admin
        :return:
        """
        return "{} {}".format(self.name, self.id)

    def active(self):
        """
        Determine whether this flowcell has been active in the 48 hours
        :return:
        """
        # time deltas are pythons measurement of time difference
        delta = datetime.timedelta(days=int(get_env_variable("MT_TIME_UNTIL_INACTIVE")))
        # If the current time minus two days is more than the last activity date, there has been no activity in 48 hours
        if (datetime.datetime.now(datetime.timezone.utc) - delta) > self.last_activity_date:

            return False
        # Activity in the last 48 hours

        return True

    class Meta:
        permissions = (
            ('view_data', 'View data'),
            ('run_analysis', 'Run analysis'),
        )


class Run(models.Model):
    """
    :purpose: Store information about each Run generated during the sequencing process in a flowcell. A flowcell can contain multiple runs, usually alternating short ones (muxscan) with long ones (real sequence data).

    Fields:

    :minion: the MinION that contains the flowcell used the produce the run data
    :flowcell: the flowcell used to produce the run data
    :owner: the run's owner; usually provides the token used in the minotour client
    :name: defines the run's name; this information usually come from the -k parameter in the minotour client
    :runid: unique identified generated by MinKnow
    :is_barcoded: true if the run has barcoded read, false otherwise
    :has_fastq: true if the sequence and quality were sent by the minotour client (stored at FastqReadExtra)
    :active: true if the Run is active; defined by a set of criteria and updated by a task
    :to_delete: to be removed
    :start_time: returns the time of the earliest read; updated by web.tasks.update_run_start_time
    """

    minion = models.ForeignKey(

        Minion,
        blank=True,
        null=True,
        related_name='runs',
        on_delete=models.CASCADE,
    )

    flowcell = models.ForeignKey(

        Flowcell,
        related_name='runs',
        on_delete=models.CASCADE,

    )

    owner = models.ForeignKey(

        settings.AUTH_USER_MODEL,
        related_name='runs',
        on_delete=models.DO_NOTHING,

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
                if int(max_channel) <= 126:
                    return 126
                elif 126 < int(max_channel) <= 512:
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


class MinionControl(models.Model):
    """
    :purpose: Store a new job for the minion to pickup and perform, example stop_minion to stop run

    Fields:
    :minion: The minion name as stored in the Minion table
    :owner: The Owner as the user
    :job: The Job that you want the minion to perform
    :custom: The contents of a custom message
    :complete: Whether the job is complete or not
    :created_date: The date that the job was created
    :modified_date: The date that this as last modified
    """
    minION = models.ForeignKey(Minion, blank=True, null=True, related_name='minioncontrol',
                               on_delete=models.CASCADE,
                               )
    owner = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='controlcontrol',
                              on_delete=models.DO_NOTHING,
                              )
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


class MinionInfo(models.Model):
    """
    :purpose: To store data on the current status of a specific MinION at a given moment in time. Doesn't matter about whether or not we have a run on.
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
    :minKNOW_warnings: If minKNOW is about to shutdown a warning is added here.
    :minknow_version: Version of minknow
    """
    minion = models.OneToOneField(
        Minion,
        related_name='currentdetails',
        on_delete=models.CASCADE,

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
    minknow_version = models.CharField(
        max_length=64,
        null=True
    )
    last_modified = models.DateTimeField(
        auto_now=True
    )

    class Meta:
        verbose_name = 'MinION Info'
        verbose_name_plural = 'MinION Info'

    def __str__(self):
        return "{} {}".format(self.minion, self.minKNOW_status)


class MinionRunStats(models.Model):
    """
    Model to store all the statistics about flowcell and MinKNOW that change during a run. Monitored once a minute,
    provides historical record of the flowcell.

    Fields
    ------

    minion: Foreign key to link to a specific MinION run.
    run_id: Foreign key to link to a specific run id.
    sample_time: Date Time field identifying the moment of data collection.
    event_yield: Sequencing yield as approximated in events by MinKNOW.
    asic_temp: Current temperature of the ASIC.
    heat_sink_temp: Current temperature of the heat sink.
    voltage_value: Current set voltage on the device.
    mean_ratio: A ratio of currents for in strand versus open pore.
    open_pore: Count of pore classification type (as named)
    in_strand: Count of pore classification type (as named)
    multiple: Count of pore classification type (as named)
    unavailable: Count of pore classification type (as named)
    unknown: Count of pore classification type (as named)
    adapater: Count of pore classification type (as named)
    pending_mux_change: Count of pore classification type (as named)
    unclassified: Count of pore classification type (as named)
    below: Count of pore classification type (as named)
    unblocking: Count of pore classification type (as named)
    above: Count of pore classification type (as named)
    good_single: Count of pore classification type (as named)
    saturated: Count of pore classification type (as named)
    inrange: Count of pore classification type (as named)
    strand: Count of pore classification type (as named)
    minKNOW_read_count: Count of reads seen at this point.
    minKNOW_histogram_values: Reporting the values of the histogram from MinKNOW
    minKNOW_histogram_bin_width: Measure of the bin width for the above histogram.
    created_date: Created Date.
    actual_max_val: Maximum read length
    n50_data: Estimated N50
    estimated_selected_bases: int
    basecalled_bases
    basecalled_fail_read_count
    basecalled_pass_read_count

    """

    minion = models.ForeignKey(

        Minion,
        related_name='currentrunstats',
        on_delete=models.CASCADE,
        null=True


    )

    run = models.ForeignKey(

        Run,
        related_name='RunStats',
        on_delete=models.CASCADE,

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

    no_pore = models.IntegerField(

        default=0
    )

    zero = models.IntegerField(

        default=0
    )

    pore = models.IntegerField(

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

    actual_max_val = models.IntegerField(
        default=0,
        null=True
    )

    n50_data = models.IntegerField(
        default=0,
        null=True
    )

    estimated_selected_bases = models.BigIntegerField(
        default=0,
        null=True
    )

    basecalled_bases = models.BigIntegerField(
        default=0,
        null=True
    )

    basecalled_fail_read_count = models.IntegerField(
        default=0,
        null=True
    )

    basecalled_pass_read_count = models.IntegerField(
        default=0,
        null=True
    )

    class Meta:
        verbose_name = 'MinION Run Stats'
        verbose_name_plural = 'MinION Run Stats'

    def __str__(self):
        return "{} {} {}".format(self.minion, self.run_id, self.sample_time)

    ## This is something to look at for optimisation
    def occupancy(self):
        if self.strand > 0 and (self.inrange > 0 or self.pore > 0):
            occupancy = round(((self.strand + self.adapter) / (self.strand + self.adapter + self.good_single + self.pore)) * 100)
        else:
            occupancy = 0
        return occupancy

    def has_basecalling(self):
        """
        Return if the run has basecalling enabled, so we display those values on the live event page
        Returns
        -------
        bool
            Whether basecalling is switched on or off
        """
        return self.basecalled_bases > 0


class MinionRunInfo(models.Model):
    """
    Model to store the data for a Run if we are monitoring MinKNOW using minFQ. Created once on Run creation. Contains all
    metadata that is fixed during a run.

    Fields
    ------
    minION
    minKNOW_current_script
    minKNOW_sample_name
    minKNOW_exp_script_purpose
    minKNOW_flow_cell_id
    minKNOW_version
    minKNOW_run_name
    run_id
    minKNOW_hash_run_id
    minKNOW_script_run_id
    minKNOW_real_sample_rate
    minKNOW_asic_id
    minKNOW_start_time
    minKNOW_colours_string
    minKNOW_computer
    experiment_type
    experiment_id
    TODO UNTIL WELLS PER CHANNEL NONE OF THESE ARE UPLOADED
    fast5_output_fastq_in_hdf
    fast5_raw
    fast5_reads_per_folder
    fastq_enabled
    fastq_reads_per_file
    filename
    flowcell_type
    kit_classification
    local_basecalling
    sample_frequency
    sequencing_kit
    user_filename_input
    wells_per_channel
    target_temp: Target temperature of minion device
    """
    minion = models.ForeignKey(

        Minion,
        related_name='currentrundetails',
        on_delete=models.CASCADE,

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

    run = models.ForeignKey(

        Run,
        related_name='RunDetails',
        on_delete=models.CASCADE,

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

    wells_per_channel = models.IntegerField(
        default=4,
        null=True
    )

    read_length_type = models.CharField(
        default="ESTIMATED BASES",
        max_length=32,
        null=True
    )

    target_temp = models.IntegerField(
        default=35,
        null=-True
    )
    class Meta:
        unique_together = ("minion", "minKNOW_hash_run_id", "minKNOW_exp_script_purpose")
        verbose_name = 'MinION Run Information'
        verbose_name_plural = 'MinION Run Information'

    def __str__(self):
        return "{} {} {}".format(self.minion, self.minKNOW_current_script, self.run_id)

    def minION_name(self):
        return self.minion.minION_name


class MinionEventType(models.Model):
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


class MinionEvent(models.Model):
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

        Minion, related_name='events',
        on_delete=models.CASCADE,

    )

    event = models.ForeignKey(

        MinionEventType,
        on_delete=models.CASCADE,
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


class MinionScripts(models.Model):
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

        Minion,
        related_name='scripts',
        on_delete=models.CASCADE,
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

    base_calling = models.BooleanField(

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


class MinionMessage(models.Model):
    """
    :purpose: Store all the Messages that the minion sends over

    Fields:
    :minion: FK linking to the minion details stored in the devices_minion table
    :run: FK linking to the Run detail
    :message: The message things like run started etc.
    :identifier: TODO matt
    :severity: Severity level of the message, provided by minKnow. 1 - normal , 2 - warning, 3 - error
    :timestamp:  The timestamp of whenn the message was sent
    """
    minion = models.ForeignKey(
        Minion,
        related_name='messages',
        on_delete=models.CASCADE,
    )

    run = models.ForeignKey(
        Run,
        related_name='runmessages',
        blank=True,
        null=True,
        on_delete=models.CASCADE,
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

    full_text = models.TextField(
        blank=True,
        default=""
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
