FROM ubuntu:latest
MAINTAINER Roberto Santos

ENV PYTHONUNBUFFERED 1

RUN mkdir -p /var/lib/minotour/apps/minotourapp /var/lib/minotour/logs /var/lib/minotour/data

WORKDIR /var/lib/minotour/apps/minotourapp

COPY . /var/lib/minotour/apps/minotourapp/

RUN apt-get update && apt-get install -y git build-essential libz-dev python3 python3-pip libmysqlclient-dev

RUN pip3 install -r requirements.txt

RUN cd extra/centrifuge-1.0.4-beta && make && cd ..

RUN cd extra/miniasm-0.3 && make && cd ..

RUN chmod -R 755 extra/*

EXPOSE 8000

CMD ["sh", "entrypoint-celery-worker.sh"]
