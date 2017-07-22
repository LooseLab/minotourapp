from celery import task
from django.core.mail import send_mail
from reads.models import MinIONmessages

@task
def contact_user(user):
    """
    At this point we would contact a user if we knew what we were doing.

    :param user:
    :return:
    """
    print (user)
    return True