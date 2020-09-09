from django.apps import AppConfig


class ReadsConfig(AppConfig):
    name = 'reads'
    def ready(self):
        import reads.services

