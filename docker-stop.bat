#!/bin/sh
# This script stops the containers used by minotour.
# Author: Roberto Santos


docker stop celery-worker-minotour-instance;

docker stop celery-minotour-instance;

docker stop web-minotour-instance;

docker stop db-minotour-instance;

docker stop redis-minotour-instance;
