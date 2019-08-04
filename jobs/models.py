from django.db import models

from reads.models import Run, Flowcell
from reference.models import ReferenceInfo


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

    run = models.ForeignKey(

        Run,
        on_delete=models.DO_NOTHING,
        related_name='runjobs',
        null=True,
        blank=True,
    )

    flowcell = models.ForeignKey(

        Flowcell,
        on_delete=models.CASCADE,
        related_name='flowcelljobs',
        null=True,
        blank=True,
    )

    job_type = models.ForeignKey(

        JobType,
        on_delete=models.DO_NOTHING,
        related_name='taskname',
    )

    reference = models.ForeignKey(

        ReferenceInfo,
        on_delete=models.CASCADE,
        related_name='referencejob',
        null=True,
        blank=True,
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

            return "{} {} {} {}".format(self.run, self.job_type, self.run.id, self.id)


        return "{} {} {}".format(self.flowcell, self.job_type, self.flowcell.id)

