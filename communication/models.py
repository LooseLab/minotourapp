from django.db import models

from django.conf import settings


class Message(models.Model):

    recipient = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        related_name='messages_recipient',
    )

    sender = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        related_name='messages_sender',
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

    def __str__(self):
        return self.title

