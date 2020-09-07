"""
Services to help run some of the tasks in the reads app.
"""
import datetime
from textwrap import wrap

from celery.utils.log import get_task_logger
from django.db.models.signals import post_save
from django.dispatch import receiver

from communication.models import NotificationConditions
from reads.models import (
    Flowcell,
    Run,
    MinionMessage,
)
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
        instance: reads.models.Flowcell
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
            **{#
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
    instance: reads.models.Run
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


@receiver(post_save, sender=MinionMessage)
def new_minion_message(sender, instance=None, created=False, **kwargs):
    """
    Save a message to the User if the minion message that has come in is a Mux or space warning.

    Parameters
    ----------
    sender:
        Unknown
    instance: reads.models.MinionMessage
        The instance of the MinionMessage that is being saved.
    created: bool
        Whether the save that is being called is creating the object in the database.
    kwargs: dict
        Key word arguments, if any.

    Returns
    -------
    None

    """
    if created:
        user = instance.minion.owner
        flowcell = instance.run.flowcell
        queryset = NotificationConditions.objects.filter(
            flowcell=flowcell, completed=False
        )
        # Messages sent as Warnings (Severity 2), Messages sent as Errors (Severity 3)
        if int(instance.severity) > 1:
            queryset = queryset.filter(notification_type="waer")
        # Mux messages me thinks
        elif instance.message.startswith("Flow cell"):
            queryset = queryset.filter(notification_type="mux")
        # Messages sent by minotTour
        elif instance.message.startswith("minoTour"):
            queryset = queryset.filter(notification_type="mino")

        if queryset:
            title = "{} from computer {} at {}".format(
                instance.message, instance.minion.computer(), datetime.datetime.now(),
            )
            chunks = wrap(title, 512)
            for chunk in chunks:
                create_and_save_message(
                    **{
                        "created": True,
                        "user": user,
                        "title": chunk,
                        "run": instance.run,
                        "flowcell": flowcell,
                    }
                )
