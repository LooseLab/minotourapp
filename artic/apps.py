from django.apps import AppConfig


class ArticConfig(AppConfig):
    """
    The config file for the Artic app.
    """
    name = 'artic'

    def ready(self):
        """
        Hook that is called when the Artic app is initialised and ready. Code below is ready.
        Returns
        -------

        """
        from artic.utils import check_artic_static_exists
        check_artic_static_exists()
