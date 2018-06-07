FROM python:latest
MAINTAINER Roberto Santos

ENV PYTHONUNBUFFERED 1
ENV MT_SECRET_KEY=''
ENV MT_DJANGO_DEBUG=''
ENV MT_DB_USER='dba'
ENV MT_DB_PASS='dba'
ENV MT_DB_HOST='mysql'
ENV MT_DB_PORT='3306'
ENV MT_DB_NAME='minotour'
ENV MT_MAILGUN_ACCESS_KEY=''
ENV MT_MAILGUN_SERVER_NAME=''
ENV MT_TWITTOKEN=''
ENV MT_TWITTOKEN_SECRET=''
ENV MT_TWITCONSUMER_KEY=''
ENV MT_TWITCONSUMER_SECRET=''
ENV MT_REFERENCE_LOCATION='/home/mbzros/a/data'
ENV PATH="$HOME/a/opt/minimap2:$HOME/a/opt/bwa:$PATH"
ENV MT_MINIMAP2="~/a/opt/minimap2/minimap2"
ENV MT_BWA="~/a/opt/bwa/bwa"
ENV MT_DB_ENGINE="django.db.backends.mysql"
#ENV MT_DB_ENGINE="django.db.backends.postgresql"

RUN mkdir -p /opt/minotour
COPY . /opt/minotourapp
WORKDIR /opt/minotourapp
RUN pip install -r requirements.txt
EXPOSE 8000
ENTRYPOINT ["python", "manage.py", "runserver", "0.0.0.0:8000"]
