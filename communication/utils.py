from communication.models import Message


def send_message(recipients, sender, title, message):

    for recipient in recipients:

        Message.objects.create(
            recipient=recipient,
            sender=sender,
            title=title,
            content=message
        )
