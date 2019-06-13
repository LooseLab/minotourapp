FROM ubuntu:latest
MAINTAINER Roberto Santos

ENV PYTHONUNBUFFERED 1
ENV MT_DB_ENGINE="django.db.backends.mysql"
ENV MT_DB_NAME="minotour"
ENV MT_DB_USER="root"
ENV MT_DB_PASS="test"
ENV MT_DB_PORT="3306"
ENV MT_DB_HOST="db-minotour-instance"
ENV MT_SECRET_KEY=''
ENV MT_DJANGO_DEBUG='True'
ENV MT_MAILGUN_ACCESS_KEY=''
ENV MT_MAILGUN_SERVER_NAME=''
ENV MT_TWITTOKEN=''
ENV MT_TWITTOKEN_SECRET=''
ENV MT_TWITCONSUMER_KEY=''
ENV MT_TWITCONSUMER_SECRET=''
ENV MT_REFERENCE_LOCATION='/var/lib/minotour/data'
ENV MT_MINIMAP2="/var/lib/minotour/apps/minotourapp/extra/minimap2/minimap2"
ENV MT_CELERY_BROKER_URL='redis://minotour_redis:6379/0'
ENV MT_CELERY_RESULT_BACKEND='redis://minotour_redis:6379/0'

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
