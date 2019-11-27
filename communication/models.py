import uuid

from django.conf import settings
from django.db import models

from reads.models import Flowcell


class Message(models.Model):

    uuid = models.UUIDField(
        default=uuid.uuid4,
        editable=False
    )

    recipient = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        related_name='messages_recipient',
        on_delete=models.DO_NOTHING,
    )

    sender = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        related_name='messages_sender',
        null=True,
        blank=True,
        on_delete=models.DO_NOTHING,
    )

    title = models.CharField(
        max_length=128,
    )

    content = models.TextField(

    )

    created_date = models.DateTimeField(
        auto_now_add=True
    )

    delivered_date = models.DateTimeField(
        null=True,
        blank=True,
    )

    read_date = models.DateTimeField(
        null=True,
        blank=True,
    )

    is_read = models.BooleanField(
        default = False,
    )

    def __str__(self):
        return self.title

    def uuidstr(self):
        return str(self.uuid)


class NotificationConditions(models.Model):
    """
        The conditions that the User has chosen to be notified about via twitter once fulfilled.
    """

    NOTIFICATION_CHOICES = [("cov", "Coverage"), ("mux", "Mux Scan Results"), ("targ", "New Target Added")]

    notification_type = models.CharField(
        max_length=4,
        choices=NOTIFICATION_CHOICES
    )

    flowcell = models.ForeignKey(
        Flowcell,
        on_delete=models.CASCADE,
        related_name="flowcell_for_conditon"
    )

    creating_user = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.CASCADE,
        related_name="condition_user"
    )

    completed = models.BooleanField(
        default=False
    )

    repeat = models.BooleanField(
        default=True
    )

    last_minKnow_message_id = models.IntegerField(
        blank=True,
        null=True
    )

    coverage = models.IntegerField(
        blank=True,
        null=True
    )
