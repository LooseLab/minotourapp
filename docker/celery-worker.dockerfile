FROM ubuntu:latest
MAINTAINER Roberto Santos

ENV PYTHONUNBUFFERED 1

RUN mkdir -p /var/lib/minotour/apps/minotourapp /var/lib/minotour/logs /var/lib/minotour/data

WORKDIR /var/lib/minotour/apps/minotourapp

COPY requirements.txt /var/lib/minotour/apps/minotourapp/

RUN apt-get update && apt-get install -y git build-essential libz-dev python3 python3-pip mariadb-client-10.1 libmysqlclient-dev

RUN pip3 install -r /var/lib/minotour/apps/minotourapp/requirements.txt

COPY docker/entrypoint.sh /var/lib/minotour/apps/minotourapp/docker/

COPY . /var/lib/minotour/apps/minotourapp/

RUN cd extra/centrifuge-1.0.4-beta && make && cd ..

RUN cd extra/miniasm-0.3 && make && cd ..

RUN chmod -R 755 extra/*

CMD ["celery", "-A", "minotourapp", "worker", "-l", "info"]