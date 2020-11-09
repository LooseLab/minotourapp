FROM python:3.8
MAINTAINER Rory Munro

ENV PYTHONUNBUFFERED 1

RUN mkdir -p /var/lib/minotour/apps/minotourapp /var/lib/minotour/logs /var/lib/minotour/data

WORKDIR /var/lib/minotour/apps/minotourapp

COPY requirements.txt /var/lib/minotour/apps/minotourapp/

RUN apt-get update && apt-get install -y python3-pip mariadb-client

RUN pip3 install -r /var/lib/minotour/apps/minotourapp/requirements.txt

COPY . /var/lib/minotour/apps/minotourapp/
