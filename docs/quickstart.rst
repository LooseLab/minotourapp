#########################
Quickstart use cheatsheet
#########################

------
Docker
------

If the docker containers have been setup following the instructions here - :doc:`docker`

Then to start the containers after restarting a computer, on unix simply run::
    sudo ./docker-start.sh

And on windows simply run::
    ./docker-start.bat

To stop containers - run the respective docker-stop file.

-----------
Development
-----------

If the minotourapp has been installed correctly, simply run::

    python3 manage.py runserver 8100

And in a seperate terminal::

    cd /path/to/minotour/code/ && source minotourenv/bin/activate && celery -A minotourapp worker -l info -B
