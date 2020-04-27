"""Apps config file
"""
from django.apps import AppConfig


class WebConfig(AppConfig):
    """
    The config file for the web app.
    """
    name = 'web'

    def ready(self):
        """
        Hook that is called when the Web app is initialised and ready. Code below is ready.
        Returns
        -------

        """
        from web.utils import validate_envs
        validate_envs()
