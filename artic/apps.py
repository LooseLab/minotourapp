from pathlib import Path

import git
from django.apps import AppConfig
from git import Repo


class ArticConfig(AppConfig):
    """
    The config file for the Artic app.
    """

    name = "artic"

    def ready(self):
        """
        Hook that is called when the Artic app is initialised and ready. Code below is ready.
        Returns
        -------

        """
        from artic.utils import check_artic_static_exists

        check_artic_static_exists()
        from minotourapp.utils import get_env_variable

        MT_VoC_PATH = get_env_variable("MT_VoC_PATH")
        if Path(f"{MT_VoC_PATH}").exists():
            print("VoC Path Found")
            ##Check if
            # cloned_repo = Repo.clone(os.path.join("https://github.com/phe-genomics/variant_definitions", Path(f"{MT_VoC_PATH}")))
            if Path(f"{MT_VoC_PATH}/variant_definitions/").exists():
                # already cloned so....
                print("Updating path")
                try:
                    repo = Repo(Path(f"{MT_VoC_PATH}/variant_definitions/"))
                    print(repo.remotes.origin.pull())
                    pass
                except git.GitCommandError as e:
                    print(
                        f"Git error, presumably being updated simultaneously {repr(e)}"
                    )
            else:
                try:
                    cloned_repo = Repo.clone_from(
                        "https://github.com/phe-genomics/variant_definitions",
                        f"{MT_VoC_PATH}/variant_definitions/",
                    )
                except git.GitCommandError as e:
                    print(
                        f"Git error, presumably being updated simultaneously {repr(e)}"
                    )
