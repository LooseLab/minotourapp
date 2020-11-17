FROM ubuntu:latest
MAINTAINER Rory Munro

ENV PYTHONUNBUFFERED 1

RUN mkdir -p /var/lib/minotour/apps/minotourapp /var/lib/minotour/logs /var/lib/minotour/data

WORKDIR /var/lib/minotour/apps/minotourapp

COPY requirements.txt /var/lib/minotour/apps/minotourapp/

RUN apt-get update && apt-get install -y git build-essential libz-dev python3 libssl-dev python3-pip mariadb-client libmysqlclient-dev

RUN pip3 install -r /var/lib/minotour/apps/minotourapp/requirements.txt

COPY docker/entrypoint.sh /var/lib/minotour/apps/minotourapp/docker/

COPY . /var/lib/minotour/apps/minotourapp/

RUN chmod -R 775 extra/*
