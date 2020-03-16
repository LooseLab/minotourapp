"""
Task_send_message.py. Code for sending tweets once a certain condition has been met on the flowcell.
"""
import pandas as pd
from celery import task
from twitter import Twitter, OAuth

from alignment.models import PafSummaryCov
from communication.models import NotificationConditions, Message
import datetime
from django.conf import settings
from celery.utils.log import get_task_logger
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


def check_coverage(flowcell, target_coverage, reference_line_id, barcode):
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

    queryset = PafSummaryCov.objects.filter(job_master__flowcell=flowcell)

    if not queryset:
        logger.warning(
            f"No PafSummaryCov objects found for flowcell with id {flowcell.id}."
        )
        return

    df = pd.DataFrame.from_records(
        queryset.values(
            "barcode_name",
            "chromosome__line_name",
            "total_length",
            "reference_line_length",
            "chromosome_id",
        )
    )

    df["coverage"] = (
        df["total_length"].div(df["reference_line_length"]).round(decimals=3)
    )
    logger.info(df.head())
    df = df.drop(columns=["reference_line_length", "total_length"])

    ref_df = df[df["chromosome_id"] == reference_line_id]
    logger.info((flowcell, target_coverage, reference_line_id, barcode))
    logger.info(ref_df)

    # ref_df = ref_df[ref_df["barcode_name"] == barcode]

    logger.info(ref_df)

    ref_df["coverage_reached"] = ref_df["coverage"] > int(target_coverage)

    results = ref_df.to_dict(orient="records")

    return results


@task()
def send_messages():
    """
    Send any unsent Messages stored in the database to the destined user using the Twitter API
    Returns
    -------
    none

    """
    time_limit = datetime.datetime.now(datetime.timezone.utc) - datetime.timedelta(hours=1)
    
    new_messages = Message.objects.filter(delivered_date=None)

    for new_message in new_messages:
        # print('Sending message: {}'.format(new_message))

        message_sent = False

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


            # status = '@{} {}'.format(new_message.recipient.extendedopts.twitterhandle,new_message.title)
            # t.statuses.update(
            #    status=status
            # )

            message_sent = True

            time.sleep(1)

        if message_sent:
            new_message.delivered_date = datetime.datetime.now(tz=datetime.timezone.utc)
            new_message.save()


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
        # Users twitter handle
        print(condition.notification_type)

        twitter_details = UserOptions.objects.get(owner=condition.creating_user)

        if not twitter_details:
            logger.info(
                f"User {condition.creating_user} does not have twitter details."
            )

        twitter_permission = twitter_details.tweet

        if not not not twitter_permission:
            logger.warning(f"This user does not have twitter permissions granted.")
            continue

        if condition.notification_type == "cov":
            # Coverage condition
            logger.info("S'up in the cov hood")
            if not condition.coverage_target:
                logger.warning(
                    f"Condition {condition.id} of type {condition.notification_type}"
                    f" should have a target coverage. Please investigate. Skipping... "
                )
                continue

            chromosome = condition.chromosome

            reference = condition.conditional_reference

            flowcell = condition.flowcell

            barcode_name = condition.barcode.name

            coverage_reached = check_coverage(
                flowcell, int(condition.coverage_target), chromosome.id, barcode_name
            )

            logger.debug("coverage reached!")
            logger.debug(coverage_reached)

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
                    print(message_text)

                    message = Message(
                        recipient=condition.creating_user,
                        sender=condition.creating_user,
                        title=message_text,
                        flowcell=flowcell,

                    )

                    message.save()

                    condition.completed = True

                    condition.save()

        if condition.notification_type in ["occu", "volt"]:
            # Count for indexing as no negative indexing
            num_entries = MinionRunStats.objects.filter(
                run_id__flowcell=condition.flowcell
            ).count()
            # Limit for condition
            upper_limit = condition.upper_limit
            lower_limit = condition.lower_limit

            if condition.notification_type == "volt":
                queryset = MinionRunStats.objects.filter(
                    run_id__flowcell=condition.flowcell
                ).order_by("-id")[:10].values_list(
                    "voltage_value", flat=True
                )
            else:
                # Else occupancy
                queryset = MinionRunStats.objects.filter(
                    run_id__flowcell=condition.flowcell
                ).order_by("-id")[:10]
                occupancy_list = []
                # Get the occupancy
                for run_stat in queryset:
                    occupancy_list.append(run_stat.occupancy())

                queryset = occupancy_list

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
                    flowcell=condition.flowcell
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
                    flowcell=condition.flowcell
                )
                message.save()
