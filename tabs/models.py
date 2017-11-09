from django.db import models
from reads.models import MinIONRun


# Create your models here.
class RunTabs(models.Model):
    class Meta:
        verbose_name_plural = 'tabs'
        # unique_together = ('tabs', 'name')

    name = models.CharField(max_length=64)

    def __str__(self):
        return self.name


class RunID(models.Model):
    class Meta:
        verbose_name_plural = 'Runs'

    run_id = models.ForeignKey(MinIONRun, on_delete=models.CASCADE, related_name='tabrunid')
    tab_id = models.ForeignKey(RunTabs, on_delete=models.CASCADE, related_name='tabs')

    def __str__(self):
        return "{}".format(self.run_id)
