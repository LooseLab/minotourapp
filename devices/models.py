from django.conf import settings
from django.db import models


class Flowcell(models.Model):

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

    minION_name = models.CharField(
        max_length=64
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
