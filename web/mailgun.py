from django.core.mail.backends.base import BaseEmailBackend


class MailgunBackend(BaseEmailBackend):
    def close(self):
        super().close()

    def open(self):
        super().open()

    def send_messages(self, email_messages):
        pass

