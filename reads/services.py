"""
Services to help run some of the tasks in the reads app.
"""
import datetime

from celery.utils.log import get_task_logger
from django.db.models.signals import post_save
from django.dispatch import receiver

from minknow_data.models import Run, Flowcell
from reads.utils import create_and_save_message

logger = get_task_logger(__name__)


@receiver(post_save, sender=Flowcell)
def new_flowcell_message(sender, instance=None, created=False, **kwargs):
    """
        If a new Flowcell is created in the database
        Parameters
        ----------
        sender:
            Unknown
        instance: minknow_data.models.Flowcell
            The instance of the flowcell that is being saved.
        created: bool
            Whether the save that is being called is creating the object in the database.
        kwargs: dict
            Key word arguments, if any.

        Returns
        -------
        None

    """
    if created:
        create_and_save_message(
            **{
                "created": True,
                "run": None,
                "user": instance.owner,
                "title": f"New flowcell {instance.name} created with sample name "
                f"{instance.sample_name} at {datetime.datetime.now()}",
                "flowcell": instance,
            }
        )


@receiver(post_save, sender=Run)
def new_run_message(sender, instance=None, created=False, **kwargs):
    """
    If a new Run is created in the database
    Parameters
    ----------
    sender:
        Unknown
    instance: minknow_data.models.Run
        The instance of the Run that is being saved.
    created: bool
        Whether the save that is being called is creating the object in the database.
    kwargs: dict
        Key word arguments, if any.

    Returns
    -------
    None

    """
    if created:
        create_and_save_message(
            **{
                "created": True,
                "user": instance.owner,
                "title": f"New run {instance.name} created on flowcell "
                f"{instance.flowcell.name} at {datetime.datetime.now()}.",
                "run": instance,
                "flowcell": instance.flowcell,
            }
        )


