FROM python:3.6
MAINTAINER Roberto Santos

ENV PYTHONUNBUFFERED 1

RUN mkdir -p /var/lib/minotour/apps/minotourapp /var/lib/minotour/logs /var/lib/minotour/data

WORKDIR /var/lib/minotour/apps/minotourapp

# COPY requirements.txt /var/lib/minotour/apps/minotourapp/

# RUN apt-get update && apt-get install -y python3-pip libmysqlclient-dev

# RUN pip3 install -r /var/lib/minotour/apps/minotourapp/requirements.txt

# COPY entrypoint.sh /var/lib/minotour/apps/minotourapp/

COPY . /var/lib/minotour/apps/minotourapp/

#RUN python3 manage.py makemigrations
#RUN python3 manage.py migrate
#RUN python3 manage.py loaddata fixtures/auxiliary_data.json

EXPOSE 8000

#CMD ["python3", "manage.py", "runserver", "0.0.0.0:8000"]
CMD ["sh", "entrypoint.sh"]
