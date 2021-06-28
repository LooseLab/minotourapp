from django.apps import AppConfig


class ReadsConfig(AppConfig):
    name = 'reads'
    def ready(self):
        # if docker instances shut with harvesting 1 we need to set it to 0
        from minotourapp.redis import redis_instance
        redis_instance.set("harvesting", 0)

