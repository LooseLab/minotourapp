FROM python:3.6
MAINTAINER Roberto Santos

ENV PYTHONUNBUFFERED 1

RUN mkdir -p /var/lib/minotour/apps/minotourapp /var/lib/minotour/logs /var/lib/minotour/data

WORKDIR /var/lib/minotour/apps/minotourapp

COPY requirements.txt /var/lib/minotour/apps/minotourapp/

RUN apt-get update && apt-get install -y python3-pip mariadb-client-10.1

RUN pip3 install -r /var/lib/minotour/apps/minotourapp/requirements.txt

COPY ./docker/entrypoint.sh /var/lib/minotour/apps/minotourapp/docker/

COPY . /var/lib/minotour/apps/minotourapp/
