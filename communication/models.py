import uuid

from django.conf import settings
from django.db import models
from django.db.models import UniqueConstraint

from reads.models import Flowcell, Barcode, Run
from reference.models import ReferenceInfo, ReferenceLine


class Message(models.Model):
    """
    Helpful docstring for bad names for the messages
    Fields
    ------
    uuid: uuid4
        A unique identifier for this message
    recipient: django.contrib.auth.models.User
        The recipient of the tweet
    sender: django.contrib.auth.models.User
        The sending User/ so the owner?
    title: str
        Twitter DM text content somehow
    content: str
        Nothing. Presumably for emails
    created_date: datetime.datetime
        The date this message was created
    delivered_date: datetime.datetime
        Date message was sent
    read_date: datetime.datetime
        No idea. Presumably when the user read it, but don't see when this was set
    is_read: bool
        Whether read or not. No idea where this is used either
    """
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

    run = models.ForeignKey(
        Run,
        null=True,
        default=None,
        on_delete=models.CASCADE
    )

    flowcell = models.ForeignKey(
        Flowcell,
        null=True,
        default=None,
        on_delete=models.CASCADE
    )

    title = models.CharField(
        max_length=512,
    )

    content = models.TextField(
        null=True,
        blank=True,
        default=""
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

        Attributes
        ---------
        notification_type: str
            A choicefield of notification type
        flowcell: reads.models.Flowcell
            The flowcell that this condition is created against
        creating_user: django.contrib.auth.models.User
            The user that created this notification. Could be linked from flowcell...
        completed: bool
            Whether the condition has been completed or not. Completed conditions are not checked again.
        repeat: bool
            Whether the condition needs to be checked even after it has been met, to see if it has been met again.
        last_minKnow_message_id: int
            The last minknow message id that we have stored.
        coverage_target: int
            If the instance is of type coverage, how much coverage we are aiming for.
        upper_limit: int
            The upper limit of the occupancy or voltage conditions.
        lower_limit: int
            The lower limit of the occupancy or voltage conditions.
        date_created: datetime.datetime
            The time the notification was created.
        reference_file: reference.models.ReferenceInfo
            The reference file that we are looking for coverage against.
        chromosome: reference.models.ReferenceLine
            The contig or chromosome from the reference that we are mapping against.
        barcode: reads.models.Barcode
            The barcode that we are mapping against.
    """

    NOTIFICATION_CHOICES = [("cov", "Coverage"), ("mux", "Mux Scan Results"), ("targ", "New Target Added"),
                            ("volt", "Voltage"), ("temp", "Temperature"),
                            ("occu", "Occupancy"), ("waer", "Warnings/Errors"), ("mino", "MinoTour"), ("sped", "Speed"),
                            ("arti", "Artic"), ("suff", "Sufficient Coverage")]

    notification_type = models.CharField(
        max_length=4,
        choices=NOTIFICATION_CHOICES
    )

    flowcell = models.ForeignKey(
        Flowcell,
        on_delete=models.CASCADE,
        related_name="flowcell_for_condition"
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

    upper_limit = models.IntegerField(
        blank=True,
        null=True,
        default=None
    )

    lower_limit = models.IntegerField(
        blank=True,
        null=True,
        default=None
    )

    date_created = models.DateField(
        auto_now_add=True
    )

    reference = models.ForeignKey(
        ReferenceInfo,
        on_delete=models.SET_NULL,
        related_name="to_chosen_reference",
        null=True,
        blank=True
    )

    chromosome = models.ForeignKey(
        ReferenceLine,
        on_delete=models.CASCADE,
        related_name="conditional_chromosome",
        null=True,
        default=None,
        blank=True
    )

    barcode = models.ForeignKey(
        Barcode,
        on_delete=models.CASCADE,
        related_name="conditional_barcode",
        null=True,
        default=None
    )

    last_minions_stats_id = models.BigIntegerField(
        default=0
    )

    UniqueConstraint(
        name='unique_communication',
        fields=['chromosome', 'reference', 'notification_type', 'coverage_target'],
    )

    def __str__(self):
        """

        Returns
        -------
        str
            The string representation of the object
        """
        return f"{self.flowcell} {self.notification_type}"
