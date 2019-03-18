#!/bin/sh
# This script creates the network and volumes used by the minotour containers.
# Author: Roberto Santos

docker network create --driver bridge minotour-network

docker volume create minotour_db_datadir

# There is no need for the builds below if the images come from docker hub.

docker build -f docker/db.dockerfile -t robertosantos/db-minotour .

docker build -f docker/web.dockerfile -t robertosantos/web-minotour .

docker build -f docker/celery.dockerfile -t robertosantos/celery-minotour .

docker build -f docker/celery-worker.dockerfile -t robertosantos/celery-worker-minotour .

docker run --name redis-minotour-instance -d --network minotour-network redis;

docker run --name db-minotour-instance -d -v minotour_db_datadir:/var/lib/mysql --network minotour-network robertosantos/db-minotour

docker run --name web-minotour-instance -p 10000:8000 -d --network minotour-network robertosantos/web-minotour

docker run --name celery-minotour-instance -d --network minotour-network robertosantos/celery-minotour

docker run --name celery-worker-minotour-instance -d --network minotour-network robertosantos/celery-worker-minotour
