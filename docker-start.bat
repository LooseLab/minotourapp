#!/bin/sh
# This script starts the containers used by minotour.
# Author: Roberto Santos

docker start redis-minotour-instance

docker start db-minotour-instance

docker start web-minotour-instance

docker start celery-minotour-instance

docker start celery-worker-minotour-instance
