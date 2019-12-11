import uuid

from django.conf import settings
from django.db import models


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
