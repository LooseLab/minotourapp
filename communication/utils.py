from communication.models import Message


def send_message(recipients, title, message, sender=None):

    for recipient in recipients:

        Message.objects.create(
            recipient=recipient,
            sender=sender,
            title=title,
            content=message
        )

def message(recipients, title, message):
    for recipient in recipients:

        Message.objects.create(
            recipient=recipient,
        #    sender=sender,
            title=title,
            content=message
        )
