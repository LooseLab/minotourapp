"""
Post save services from minknow message creation. Creates a message to be tweeted
"""
from textwrap import wrap

from django.db.models.signals import post_save
from django.dispatch import receiver

from communication.models import NotificationConditions
from minknow_data.models import MinionMessage
from reads.utils import create_and_save_message


@receiver(post_save, sender=MinionMessage)
def new_minion_message(sender, instance=None, created=False, **kwargs):
    """
    Save a message to the User if the minion message that has come in is a Mux or space warning.

    Parameters
    ----------
    sender:
        Unknown
    instance: minknow_data.models.MinionMessage
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
