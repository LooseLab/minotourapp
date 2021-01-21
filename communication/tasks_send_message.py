"""
Task_send_message.py. Code for sending tweets once a certain condition has been met on the flowcell.
"""
import datetime
import time

from celery.utils.log import get_task_logger
from django.conf import settings
from django.db.models import F, Sum, Avg
from twitter import Twitter, OAuth, TwitterHTTPError, TwitterError

from alignment.models import PafSummaryCov
from artic.models import ArticBarcodeMetadata
from communication.models import NotificationConditions, Message
from minknow_data.models import MinionRunStats, MinionControl
from minotourapp.celery import app
from reads.models import UserOptions

logger = get_task_logger(__name__)


def return_tweet(
    tweet_type,
    target_coverage=None,
    reference_file_name=None,
    chromosome=None,
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
    chromosome: reference.models.ReferenceLine
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
    chromosome_name = chromosome if chromosome else "Complete genome"
    text_lookup = {
        "coverage": f"Target coverage: {target_coverage}x reached against reference file:"
        f" {reference_file_name}, chromosome name: {chromosome_name}. Recorded coverage is: {observed_coverage}x.",
        "new_target": f"New Target added on Read Until iteralign run: {target_added}",
        "volt_occu": f"Limit breached for {volt_occu} on Flowcell {flowcell}. Limit set was {limit}.",
    }
    tweet_text = text_lookup[tweet_type]
    return tweet_text


def check_coverage(flowcell, target_coverage, reference_line, barcodes):
    """
    Check the coverage of a mapping task.
    Parameters
    ----------
    flowcell: minknow_data.models.Flowcell
        The flowcell that sequenced the reads for this mapping task.
    target_coverage: int
        The desired amount of coverage to be reached
    reference_line: reference.models.ReferenceLine
        Could be null for genome level coverage, otherwise reference line object
    barcodes: list
        The primary key of the barcodes

    Returns
    -------
    results: list
         list of dictionaries of results showing whether or not we have reached coverage for the chosen chromosomes.
    """
    if reference_line:
        queryset = PafSummaryCov.objects.filter(
            job_master__flowcell=flowcell,
            chromosome_id=reference_line.id,
            coverage__gte=target_coverage,
            barocde_id__in=barcodes
        ).values("coverage")
        queryset = [
            {"coverage_reached": True, "coverage": coverage} for coverage in queryset
        ]
    else:
        queryset = (
            PafSummaryCov.objects.filter(job_master__flowcell=flowcell, barcode_id__in=barcodes)
            .values(data_name=F("reference_name"))
            .annotate(
                Sum("total_yield"),
                Sum("reference_line_length"),
                Avg("average_read_length"),
            )
        )
        queryset = [
            {
                "coverage": q["total_yield__sum"] / q["reference_line_length__sum"],
                "coverage_reached": True,
            }
            for q in queryset
            if q["total_yield__sum"] / q["reference_line_length__sum"] > target_coverage
        ]
    if not queryset:
        logger.warning(
            f"No PafSummaryCov objects found for flowcell with id {flowcell.id}."
        )
        return
    return queryset


@app.task
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
        delivered_date=datetime.datetime.now(datetime.timezone.utc)
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
            finally:
                new_message.delivered_date = datetime.datetime.now(
                    tz=datetime.timezone.utc
                )
                new_message.save()
                time.sleep(1)
                return


def coverage_notification(condition):
    """
    Run calculations for coverage statistics - this is pretty gahbahj
    Parameters
    ----------
    condition: communication.models.NotificationConditions
        The condition we are dealing with
    Returns
    -------

    """
    if not condition.coverage_target:
        logger.warning(
            f"Condition {condition.id} of type {condition.notification_type}"
            f" should have a target coverage. Please investigate. Skipping... "
        )
        return
    chromosome = condition.chromosome
    reference = condition.reference
    flowcell = condition.flowcell
    run = flowcell.runs.last()
    barcodes = list(condition.mapping_condition_barcodes.all().values_list("barcode_id", flat=True))
    coverage_reached = check_coverage(
        flowcell, int(condition.coverage_target), chromosome, barcodes
    )
    logger.debug("coverage reached!")
    logger.debug(coverage_reached)
    if not coverage_reached:
        return
    # Reached is a dictionary, one dict for each Paf Summary Cov pulled back in the check_coverage function
    for reached in coverage_reached:
        if reached.get("coverage_reached", False):
            # Let's create a message
            message_text = return_tweet(
                "coverage",
                target_coverage=int(condition.coverage_target),
                reference_file_name=reference.name,
                chromosome=chromosome,
                observed_coverage=reached.get("coverage", "Undefined"),
            )
            message = Message(
                recipient=condition.creating_user,
                sender=condition.creating_user,
                title=message_text,
                flowcell=flowcell,
            )
            message.save()
            if condition.run_until:
                create_stop_run(run.minion)
            condition.completed = True
            condition.save()


def create_stop_run(minion):
    """
    Create a stop run minion control command that is picked up by minFQ.
    Parameters
    ----------
    minion: minknow_data.models.Minion
        The minion that the run is currently active on
    Returns
    -------
    None
    """
    MinionControl.objects.create(
        custom="", job="stop_minion", minion=minion, complete=False
    )


def minion_statistic_check(condition):
    """
        Run calculations for minion statistics, such as speed, occupation, voltage, if 10 of the last 10 MinionStatitics
        are out of the limits
        on the notification condition, send tweet and if run until enabled stop run
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
    if MinionRunStats.objects.filter(run__flowcell=condition.flowcell).count():
        queryset = MinionRunStats.objects.filter(
            run__flowcell=condition.flowcell, id__gte=condition.last_minions_stats_id,
        ).order_by("-id")[:10]
    else:
        return
    last_mrs = queryset.last()
    if condition.notification_type == "volt":
        queryset = queryset.values_list("voltage_value", flat=True)
    else:
        queryset = [minion_statistic.occupancy() for minion_statistic in queryset]
    upper_list = list(filter(lambda x: x >= upper_limit, queryset))
    lower_list = list(filter(lambda x: x <= lower_limit, queryset))
    if len(upper_list) == 10 or len(lower_list) == 10:
        limit = upper_limit if len(upper_list) == 10 else lower_limit
        twitter_details = UserOptions.objects.get(owner=condition.creating_user)
        if twitter_details.twitterhandle and twitter_details.tweet:
            text = return_tweet(
                tweet_type="volt_occu",
                volt_occu="Voltage",
                limit=limit,
                flowcell=condition.flowcell.name,
            )
            message = Message(
                recipient=condition.creating_user,
                sender=condition.creating_user,
                title=text,
                flowcell=condition.flowcell,
            )
            message.save()
        if condition.run_until:
            # TODO if we have a run until bug look here as conditions should really be run aware?
            create_stop_run(last_mrs.minion)
        condition.last_minions_stats_id = queryset.last().id
        condition.completed = True
        condition.save()

@app.task
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
        if condition.notification_type in ["occu", "volt", "sped"]:
            minion_statistic_check(condition)
        twitter_details = UserOptions.objects.get(owner=condition.creating_user)
        if not twitter_details:
            logger.warning(
                f"User {condition.creating_user} does not have twitter details."
            )
        twitter_permission = twitter_details.tweet
        if not twitter_permission and twitter_details.twitterhandle:
            logger.warning(f"This user does not have twitter permissions granted.")
            continue
        if condition.notification_type == "suff":
            check_artic_sufficient_coverage(condition)
        if condition.notification_type == "arti":
            check_artic_has_fired(condition)
        if condition.notification_type == "cov":
            coverage_notification(condition)
        if condition.notification_type == "pred":
            check_artic_run_until(condition)


def check_artic_run_until(condition):
    """
    Check that Artic has finished Run Until
    Parameters
    ----------
    condition: communication.models.NotificationConditions
        The artic notification condition

    Returns
    -------

    """
    if hasattr(condition.flowcell.runs.last(), "summary"):
        if condition.flowcell.runs.last().summary.first_read_start_time < datetime.datetime.now(datetime.timezone.utc) - datetime.timedelta(hours=1):
            if ArticBarcodeMetadata.objects.filter(
                flowcell_id=condition.flowcell, projected_to_finish=True, has_finished=False
            ).count():
                m, created = Message.objects.get_or_create(
                    recipient=condition.creating_user,
                    sender=condition.creating_user,
                    title__startswith=f"Artic run on {condition.flowcell.name} has potentially finished,"
                                      f" further sequencing may be a waste. Consider stopping this run.",
                    flowcell=condition.flowcell,
                    defaults={
                        "title": f"Artic run on {condition.flowcell.name} has potentially finished,"
                                 f" further sequencing may be a waste. Consider stopping this run.",
                    },
                )
                condition.completed = True
                condition.save()



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
        m, created = Message.objects.get_or_create(
            recipient=condition.creating_user,
            sender=condition.creating_user,
            title__startswith=f"Artic pipeline has finished for barcode {barcode['barcode__name']}",
            flowcell=condition.flowcell,
            defaults={
                "title": f"Artic pipeline has finished for barcode {barcode['barcode__name']} at {time}"
            },
        )


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
    ).values("barcode__name", "average_coverage", "job_master_id")
    for barcode in queryset:
        time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        m, created = Message.objects.get_or_create(
            recipient=condition.creating_user,
            sender=condition.creating_user,
            title__startswith=f"Sufficient Coverage reached for Barcode {barcode['barcode__name']} in Artic task {barcode['job_master_id']}",
            flowcell=condition.flowcell,
            defaults={
                "title": f"Sufficient Coverage reached for Barcode {barcode['barcode__name']} in Artic task {barcode['job_master_id']} with an average coverage of {barcode['average_coverage']} at {time}."
            },
        )
