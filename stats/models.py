from django.db import models

# Create your models here.
from reads.models import Barcode, FastqReadType, Run


class RunStatisticBarcode(models.Model):

    run_id = models.ForeignKey(
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

    def __str__(self):
        return "{} {} {} {}".format(
            self.run_id,
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

    def __str__(self):
        return "{} {} {}".format(self.run_id, self.channel_number, self.read_count)


class HistogramSummary(models.Model):

    BIN_WIDTH = 900

    run_id = models.ForeignKey(
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

    def __str__(self):
        return "{} {} {}".format(self.run_id, self.read_type, self.bin_width)


class RunSummaryBarcode(models.Model):

    run_id = models.ForeignKey(
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
