#!/bin/sh

docker network create --driver bridge minotour-network;

docker volume create minotour_db_datadir;

docker build -f db.dockerfile -t robertosantos/db-minotour .

docker build -f web.dockerfile -t robertosantos/web-minotour .

docker build -f .\celery.dockerfile -t robertosantos/celery-minotour .

docker build -f .\celery-worker.dockerfile -t robertosantos/celery-worker-minotour .

docker run --name redis-minotour-instance -p 6379:6379 -d --network minotour-network redis;

docker run --name db-minotour-instance -p 3307:3306 -d -v minotour_db_datadir:/var/lib/mysql --network minotour-network robertosantos/db-minotour;

docker run --name web-minotour-instance -p 10000:8000 -d --network minotour-network robertosantos/web-minotour;

docker run --name celery-minotour-instance --network minotour-network robertosantos/celery-minotour

docker run --name celery-worker-minotour-instance --network minotour-network robertosantos/celery-worker-minotour
