from django.db import models

from reads.models import Run, Flowcell, JobMaster
from reference.models import ReferenceLine


class RejectedFastqRead(models.Model):

    """
    It holds information about reads rejected by the read until algorithm.

    1) This model can contain reads that don't exist in the reads.FastqRead model.
    """
    flowcell = models.ForeignKey(
        Flowcell,
        null=True,
        default=None,
        on_delete=models.CASCADE,
        related_name="rejected_read_flowcell"
    )

    run = models.ForeignKey(
        Run,
        null=True,
        default=None,
        on_delete=models.CASCADE,
    )

    read_id = models.CharField(

        max_length=96
    )

    def __str__(self):

        return "{} - {}".format(self.read_id, self.run.runid)


class AcceptedFastqRead(models.Model):

    """
    It holds information about reads rejected by the read until algorithm.

    1) This model can contain reads that don't exist in the reads.FastqRead model.
    """
    flowcell = models.ForeignKey(
        Flowcell,
        on_delete=models.CASCADE,
        related_name="accepted_read_flowcell"
    )

    run = models.ForeignKey(

        Run,
        on_delete=models.CASCADE,
    )

    read_id = models.CharField(

        max_length=96
    )

    def __str__(self):

        return "{} - {}".format(self.read_id, self.run.runid)


class ExpectedBenefitChromosomes(models.Model):
    """
        This contains the chromosomes we have in a n EB task to populate the select box on the expected Benefits page
    """
    task = models.ForeignKey(
        JobMaster,
        on_delete=models.CASCADE,
        related_name="expectedBenefitTask"
    )
    chromosome = models.ForeignKey(
        ReferenceLine,
        on_delete=models.CASCADE,
        related_name="expectedBenefitChromosome"
    )

    def __str__(self):

        return "Task Id - {}, Chromosome name - {}".format(self.task.id, self.chromosome.line_name)
