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



def return_tweet(tweet_type, target_coverage=None, reference_file_name=None, mux_pores_left=None, target_added=None):
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
    mux_pores_left: int
        Number of pores left functional after a mux scan. Default None.
    target_added: str
        The species added to the iteralign rejection list.

    Returns
    -------
    tweet_text: str
        The text to be contained in the twitter DM.

    """

    text_lookup = {"coverage": f"target coverage: {target_coverage}x reached against reference file: {reference_file_name}",
                   "mux": f"mux scan finished. Pores left: {mux_pores_left}",
                   "new_target": f"new Target added on Read Until iteralign run: {target_added}"}

    tweet_text = text_lookup[tweet_type]

    return tweet_text


def check_coverage(flowcell_id, target_coverage, reference_line):
    """
    Check the coverage of a mapping task.
    Parameters
    ----------
    flowcell_id: int
        The primary key of a flowcell mapping task.
    target_coverage: int
        The desired amount of coverage to be reached
    reference_line: str
        The name of the chromosome or contig to have the coverage checked.

    Returns
    -------
    coverage_reached: bool
        Whether the coverage has been reached.
    """

    queryset = PafSummaryCov.objects.filter(job_master__flowcell__id=flowcell_id)

    if not queryset:
        return

    df = pd.DataFrame.from_records(
        queryset.values('barcode_name', 'reference_line_name', 'total_length',
                        'reference_line_length'))

    df["coverage"] = df['total_length'].div(df['reference_line_length']).round(decimals=3)

    df = df.drop(columns=["reference_line_length", "total_length"])

    ref_df = df[df["reference_line_name"] == reference_line]

    ref_df[ref_df["barcode_name"]]

@task
def send_messages():
    """
    Send any unsent Messages stored in the database to the destined user using the Twitter API
    Returns
    -------

    """
    new_messages = Message.objects.filter(delivered_date=None)

    for new_message in new_messages:
        # print('Sending message: {}'.format(new_message))

        message_sent = False

        if new_message.recipient.extendedopts.tweet \
                and new_message.recipient.extendedopts.twitterhandle != '':
            TWITTOKEN = settings.TWITTOKEN
            TWITTOKEN_SECRET = settings.TWITTOKEN_SECRET
            TWITCONSUMER_KEY = settings.TWITCONSUMER_KEY
            TWITCONSUMER_SECRET = settings.TWITCONSUMER_SECRET

            t = Twitter(
                auth=OAuth(TWITTOKEN, TWITTOKEN_SECRET, TWITCONSUMER_KEY, TWITCONSUMER_SECRET)
            )

            t.direct_messages.new(
                user=new_message.recipient.extendedopts.twitterhandle,
                text=new_message.title
            )
            # status = '@{} {}'.format(new_message.recipient.extendedopts.twitterhandle,new_message.title)
            # t.statuses.update(
            #    status=status
            # )

            message_sent = True

        if message_sent:
            print('inside message_sent')
            new_message.delivered_date = datetime.now(tz=pytz.utc)
            new_message.save()

@task()
def send_tweet_upon_condition():
    """
    This function is executed by django beat. It checks the database for a set of conditions the user has chosen, and
    if the conditions are met, sends the user a tweet alerting them.

    Returns
    -------
    None
    """
    # Get the date 7 days ago
    seven_day_date_limit = datetime.date.today() - datetime.timedelta(days=7)
    # Get conditions created in the last 7 days and conditions that are incomplete
    active_conditions = NotificationConditions.objects.filter(completed=False, date_created__gte=seven_day_date_limit)


    for condition in active_conditions:
        print(condition.notification_type)
        if condition.notification_type is "cov":




