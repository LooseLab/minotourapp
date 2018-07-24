from django.db import models

from reads.models import Run


class RejectedFastqRead(models.Model):

    """
    It holds information about reads rejected by the read until algorithm.

    1) This model can contain reads that don't exist in the reads.FastqRead model.
    """

    run = models.ForeignKey(

        Run,
    )

    read_id = models.CharField(

        max_length=96
    )

    def __str__(self):

        return "{} - {}".format(self.read_id, self.run.runid)
