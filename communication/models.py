import uuid
from django.conf import settings
from django.db import models

from reads.models import Flowcell, Barcode
from reference.models import ReferenceInfo, ReferenceLine


# class ConditionalBarcodes(models.Model):
#     """
#     store the conditional Barcodes
#     """
#     link = models.ForeignKey(
#         Barcode,
#         on_delete=models.CASCADE,
#         related_name="chosen_barcodes"
#     )
#
#
# class ConditionalChromosomes(models.Model):
#     """
#     Store the conditionalChromosomes
#     """
#     link = models.ForeignKey(
#         ReferenceLine,
#         on_delete=models.CASCADE,
#         related_name="chosen_chromosomes"
#     )


class Message(models.Model):

    uuid = models.UUIDField(
        default=uuid.uuid4,
        editable=False
    )

    recipient = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        related_name='messages_recipient',
        on_delete=models.CASCADE,
        blank=True,
        null=True
    )

    sender = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        related_name='messages_sender',
        null=True,
        blank=True,
        on_delete=models.CASCADE,
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
        default=False,
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
        null=True,
        default=0
    )

    coverage_target = models.IntegerField(
        blank=True,
        null=True
    )

    date_created = models.DateField(
        auto_now_add=True
    )

    reference_file = models.ForeignKey(
        ReferenceInfo,
        on_delete=models.SET_NULL,
        name="conditional_reference",
        related_name="to_chosen_reference",
        null=True
    )

    chromosome = models.ForeignKey(
        ReferenceLine,
        on_delete=models.CASCADE,
        related_name="conditional_chromosome",
        null=True,
        default=None
    )

    barcode = models.ForeignKey(
        Barcode,
        on_delete=models.CASCADE,
        related_name="conditional_barcode",
        null=True,
        default=None
    )

    def __str__(self):
        """

        Returns
        -------
        str
            The string representation of the object
        """
        return f"{self.flowcell} {self.notification_type}"
