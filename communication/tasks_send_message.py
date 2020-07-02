"""
Task_send_message.py. Code for sending tweets once a certain condition has been met on the flowcell.
"""
import datetime
import time

from celery import task
from celery.utils.log import get_task_logger
from django.conf import settings
from twitter import Twitter, OAuth, TwitterHTTPError, TwitterError

from alignment.models import PafSummaryCov
from artic.models import ArticBarcodeMetadata
from communication.models import NotificationConditions, Message
from reads.models import UserOptions, MinionRunStats

logger = get_task_logger(__name__)


def return_tweet(
    tweet_type,
    target_coverage=None,
    reference_file_name=None,
    chromosome_name=None,
    observed_coverage=None,
    target_added=None,
    limit=None,
    volt_occu=None,
    flowcell=None,
    last_value=None,
):
    """
    Return the tweet to be sent for the
    Parameters
    ----------
    tweet_type: str
        The name of the condition that has been met.
    target_coverage: int
        The target coverage that has been reached by a mapping task. Default None.
    reference_file_name: str
        The name of the reference that the coverage has been reached against. Default None.
    chromosome_name: str
        The name of the chromosome we have recorded coverage against.
    observed_coverage: int
        Depth of coverage that we have seen across the genome.
    target_added: str
        The species added to the iteralign rejection list.
    limit: str
        The limit that we have passed for the Occupancy / Voltage notifications
    volt_occu: str
        Voltage or occupancy
    flowcell: str
        The flowcell name.
    last_value: str
        The last value for Occupancy or Voltage that was recorded.


    Returns
    -------
    tweet_text: str
        The text to be contained in the twitter DM.

    """

    text_lookup = {
        "coverage": f"Target coverage: {target_coverage}x reached against reference file:"
        f" {reference_file_name}, chromosome name: {chromosome_name}. Recorded coverage is: {observed_coverage}x.",
        "new_target": f"New Target added on Read Until iteralign run: {target_added}",
        "volt_occu": f"Limit breached for {volt_occu} on Flowcell {flowcell}. Limit set was {limit}.",
    }
    tweet_text = text_lookup[tweet_type]
    return tweet_text


def check_coverage(flowcell, target_coverage, reference_line_id, barcode=""):
    """
    Check the coverage of a mapping task.
    Parameters
    ----------
    flowcell: reads.models.Flowcell
        The flowcell that sequenced the reads for this mapping task.
    target_coverage: int
        The desired amount of coverage to be reached
    reference_line_id: id
        The name of the chromosome or contig to have the coverage checked.
    barcode: str
        The name of the barcode

    Returns
    -------
    results: list
         list of dictionaries of results showing whether or not we have reached coverage for the chosen chromosomes.
    """

    queryset = PafSummaryCov.objects.filter(
        job_master__flowcell=flowcell,
        chromosome_id=reference_line_id,
        coverage__gte=target_coverage,
    )
    if not queryset:
        logger.warning(
            f"No PafSummaryCov objects found for flowcell with id {flowcell.id}."
        )
        return
    return queryset.values("coverage")


@task()
def send_messages():
    """
    Send any unsent Messages stored in the database to the destined user using the Twitter API
    Returns
    -------
    none

    """
    time_limit = datetime.datetime.now(datetime.timezone.utc) - datetime.timedelta(
        hours=1
    )
    # Set the messages to delivered if they are more than hour old
    Message.objects.filter(created_date__lt=time_limit, delivered_date=None).update(
        delivered_date=datetime.datetime.now()
    )
    # Get the new messages
    new_messages = Message.objects.filter(delivered_date=None)
    for new_message in new_messages:
        if new_message.created_date < time_limit:
            continue
        message_sent = False
        # If the user has set up twitter options
        if (
            new_message.recipient.extendedopts.tweet
            and new_message.recipient.extendedopts.twitterhandle != ""
        ):
            TWITTOKEN = settings.TWITTOKEN
            TWITTOKEN_SECRET = settings.TWITTOKEN_SECRET
            TWITCONSUMER_KEY = settings.TWITCONSUMER_KEY
            TWITCONSUMER_SECRET = settings.TWITCONSUMER_SECRET
            t = Twitter(
                auth=OAuth(
                    TWITTOKEN, TWITTOKEN_SECRET, TWITCONSUMER_KEY, TWITCONSUMER_SECRET
                )
            )
            try:
                t.direct_messages.events.new(
                    _json={
                        "event": {
                            "type": "message_create",
                            "message_create": {
                                "target": {
                                    "recipient_id": t.users.show(
                                        screen_name=new_message.recipient.extendedopts.twitterhandle
                                    )["id"]
                                },
                                "message_data": {"text": new_message.title},
                            },
                        }
                    }
                )
            except (TwitterHTTPError, TwitterError) as e:
                logger.error(str(e))
                return
            finally:
                message_sent = True
                if message_sent:
                    new_message.delivered_date = datetime.datetime.now(
                        tz=datetime.timezone.utc
                    )
                    new_message.save()
                    time.sleep(1)


def coverage_notification(condition):
    """
    Run calculations for coverage statistics
    Parameters
    ----------
    condition: communication.models.NotificationConditions
        The condition we are dealing with
    Returns
    -------

    """
    logger.info("S'up in the cov hood")
    if not condition.coverage_target:
        logger.warning(
            f"Condition {condition.id} of type {condition.notification_type}"
            f" should have a target coverage. Please investigate. Skipping... "
        )
        return
    chromosome = condition.chromosome
    reference = condition.reference
    flowcell = condition.flowcell
    # barcode_name = condition.barcode.name

    coverage_reached = check_coverage(
        flowcell, int(condition.coverage_target), chromosome.id
    )

    logger.debug("coverage reached!")
    logger.debug(coverage_reached)
    if not coverage_reached:
        return

    # Reached is a dictionary, one dict for each Paf Summary Cov pulled back in the check_coverage function
    for reached in coverage_reached:
        if reached.get("coverage_reached", None):
            # Let's create a message
            message_text = return_tweet(
                "coverage",
                target_coverage=int(condition.coverage_target),
                reference_file_name=reference.name,
                chromosome_name=chromosome.line_name,
                observed_coverage=reached.get("coverage", "Undefined"),
            )
            message = Message(
                recipient=condition.creating_user,
                sender=condition.creating_user,
                title=message_text,
                flowcell=flowcell,
            )
            message.save()
            condition.completed = True
            condition.save()


def minion_statistic_check(condition):
    """
        Run calculations for minion statistics, such as speed, occupation, voltage
        Parameters
        ----------
        condition: communication.models.NotificationConditions
            The condition we are dealing with
        Returns
        -------

        """
    # Limit for condition
    upper_limit = condition.upper_limit
    lower_limit = condition.lower_limit
    if (
        condition.last_minions_stats_id == 0
        and MinionRunStats.objects.filter(run__flowcell=condition.flowcell).count()
    ):
        queryset = MinionRunStats.objects.filter(run__flowcell=condition.flowcell)[:10]
        condition.last_minions_stats_id = queryset[9].id
    elif not MinionRunStats.objects.filter(run__flowcell=condition.flowcell).count():
        return
    else:
        queryset = MinionRunStats.objects.filter(
            run__flowcell=condition.flowcell, id__gte=condition.last_minions_stats_id,
        )
        condition.last_minions_stats_id = queryset.last().id
    if condition.notification_type == "volt":
        queryset.values_list("voltage_value", flat=True)
    else:
        queryset = [minion_statistic.occupancy() for minion_statistic in queryset]
    upper_list = list(filter(lambda x: x >= upper_limit, queryset))
    lower_list = list(filter(lambda x: x <= lower_limit, queryset))
    if len(upper_list) > 6:
        text = return_tweet(
            "volt_occu",
            upper_limit,
            "Voltage",
            condition.flowcell.name,
            queryset[len(queryset)],
        )
        message = Message(
            recipient=condition.creating_user,
            sender=condition.creating_user,
            title=text,
            flowcell=condition.flowcell,
        )
        message.save()
    elif len(lower_list) > 6:
        text = return_tweet(
            "volt_occu",
            lower_limit,
            "Voltage",
            condition.flowcell.name,
            queryset[len(queryset)],
        )
        message = Message(
            recipient=condition.creating_user,
            sender=condition.creating_user,
            title=text,
            flowcell=condition.flowcell,
        )
        message.save()


@task()
def check_condition_is_met():
    """
    This function is executed by django beat. It checks the database for a set of conditions the user has chosen, and
    if the conditions are met, creates a message to be sent to them.

    Returns
    -------
    None
    """
    # Get the date 7 days ago
    seven_day_date_limit = datetime.date.today() - datetime.timedelta(days=7)
    # Get conditions created in the last 7 days and conditions that are incomplete
    active_conditions = NotificationConditions.objects.filter(
        completed=False, date_created__gte=seven_day_date_limit
    )
    logger.info(f"Active conditions {active_conditions}")
    for condition in active_conditions:
        twitter_details = UserOptions.objects.get(owner=condition.creating_user)
        if not twitter_details:
            logger.info(
                f"User {condition.creating_user} does not have twitter details."
            )
        twitter_permission = twitter_details.tweet
        if not not not twitter_permission:
            logger.warning(f"This user does not have twitter permissions granted.")
            continue
        if condition.notification_type == "suff":
            check_artic_sufficient_coverage(condition)
        if condition.notification_type == "arti":
            check_artic_has_fired(condition)
        if condition.notification_type == "cov":
            coverage_notification(condition)
        if condition.notification_type in ["occu", "volt", "sped"]:
            minion_statistic_check(condition)


def check_artic_has_fired(condition):
    """
    check if a barocde in the artic task has gone through the pipelin
    Parameters
    ----------
    condition: communication.models.NotificationConditions
        The condition of this type
    Returns
    -------

    """
    queryset = ArticBarcodeMetadata.objects.filter(
        flowcell_id=condition.flowcell, has_finished=True
    ).values("barcode__name")
    for barcode in queryset:
        time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        m = Message(
            recipient=condition.creating_user,
            sender=condition.creating_user,
            title=f"Artic pipeline has finished for barcode {barcode['barcode__name']} at {time}",
            flowcell=condition.flowcell,
        )
        m.save()


def check_artic_sufficient_coverage(condition):
    """
    Check if the artic coverage requiremnts are met
    Parameters
    ----------
    condition: communication.models.NotificationConditions
        The condition of this type
    Returns
    -------

    """
    queryset = ArticBarcodeMetadata.objects.filter(
        flowcell_id=condition.flowcell, has_sufficient_coverage=True
    ).values("barcode__name", "average_coverage")
    for barcode in queryset:
        m = Message(
            recipient=condition.creating_user,
            sender=condition.creating_user,
            title=f"Sufficient Coverage reached for Barocde {barcode['barcode__name']} in Artic task,"
            f" with an average coverage of {barcode['average_coverage']}",
            flowcell=condition.flowcell,
        )
        m.save()
