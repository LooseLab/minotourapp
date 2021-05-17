"""Apps config"""
from pathlib import Path

from django.apps import AppConfig
from django.conf import settings


def check_media_roots_exist():
    """
    Check the django specified media root exists, and if not creates them
    Returns
    -------

    """
    media_root = Path(settings.MEDIA_ROOT)
    if not media_root.exists():
        media_root.mkdir()
    if not (media_root / "minimap2_indexes").exists():
        (media_root / "minimap2_indexes").mkdir()
    if not (media_root / "reference_files").exists():
        (media_root / "reference_files").mkdir()

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
