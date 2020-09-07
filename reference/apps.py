"""Apps config"""
from django.apps import AppConfig

from reference.utils import check_media_roots_exist


class ReferenceConfig(AppConfig):
    """
    Config reference App
    """
    name = 'reference'

    def ready(self):
        """
        Hook that is called when the Reference app is initialised and ready. Code below is ready.
        Returns
        -------

        """
        check_media_roots_exist()
