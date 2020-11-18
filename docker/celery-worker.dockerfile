FROM ubuntu:latest
MAINTAINER Rory Munro
SHELL [ "/bin/bash", "--login", "-c" ]
ENV MT_DB_ENGINE='django.db.backends.mysql'
ENV MT_DB_HOST='localhost'
#ENV MT_DB_NAME='django_test'
ENV MT_DB_NAME='minotour'
ENV MT_DB_USER='dba'
ENV MT_DB_PASS='test'
ENV MT_DB_PORT='3306'
ENV MT_DJANGO_DEBUG=0
ENV MT_MAILGUN_ACCESS_KEY=''
ENV MT_MAILGUN_SERVER_NAME=''
ENV MT_MINIMAP2='extra/minimap2-2.17_x64-linux/minimap2'
# Write uploaded references to this directory
ENV MT_REFERENCE_LOCATION='<path/to/references/directory/>'
ENV MT_TWITCONSUMER_KEY='<Twitter consumer key>'
ENV MT_TWITCONSUMER_SECRET='<Twitter consumer secret>'
ENV MT_TWITTOKEN='<Twitter token>'
ENV MT_TWITTOKEN_SECRET='<Twitter token secret>'
# Generate a secret key with python manage.py shell -c 'from django.core.management import utils; print(utils.get_random_secret_key())'
ENV MT_SECRET_KEY='hello_secret_key!!'
ENV MT_CENTRIFUGE="extra/centrifuge-1.0.4-beta/centrifuge"
# Read centrifuge index from this directory
ENV MT_CENTRIFUGE_INDEX="/var/lib/minotour/apps/data"
# Write Celery log file in this directory
ENV MT_LOG_FOLDER='/var/lib/minotour/logs'
ENV MT_CELERY_BROKER_URL='redis://localhost:6379/0'
ENV MT_CELERY_RESULT_BACKEND='redis://localhost:6379/0'
#ENV MT_ARTIC_RESULTS_DIR='/media/rory/Bioinformatics'
# Store the Artic results in this directory
ENV MT_ARTIC_RESULTS_DIR='/var/lib/minotour/apps/data'
 # Read the Artic scheme directories, found where the Artic NCOV piprline github repo was clones. Example  Projects/artic-ncov2019/primer_schemes
ENV MT_ARTIC_SCEHEME_DIR='<Path/To/Artic/Code/Scheme_Directories>'
# Maximum number of Artic tasks to be run at any one time
ENV MT_ARTIC_MAX_PIPELINES="1"
# Number of Celery workers to start, no more than number of CPU cores
ENV MT_CELERY_WORKER_COUNT="2"
 # Number of celery minimap2 threads to start, no more than num CPU cores
ENV MT_CELERY_MINIMAP2_WORKER_COUNT="1"
 # Number of days since last flowcell activity until a flowcell is listed as inactive
ENV MT_TIME_UNTIL_INACTIVE="7"
 # Number of days since last flowcell activity until a flowcell is archived and fastq data deleted
ENV MT_TIME_UNTIL_ARCHIVE="14"
ENV MT_SKIP_SAVING_SEQUENCE=0
ARG PATH="/root/miniconda3/bin:${PATH}"
ENV PATH="/root/miniconda3/bin:${PATH}"
ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV PYTHONUNBUFFERED 1

RUN mkdir -p /var/lib/minotour/apps/minotourapp /var/lib/minotour/logs /var/lib/minotour/data /var/lib/minotour/apps/data/references /var/lib/minotour/apps/minotour_static

WORKDIR /var/lib/minotour/apps/minotourapp

COPY requirements.txt /var/lib/minotour/apps/minotourapp/

RUN apt-get update && apt-get install -y git build-essential libz-dev python3-dev libssl-dev python3-pip mariadb-client default-libmysqlclient-dev

RUN apt-get install wget
RUN wget --quiet \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \

RUN conda init bash
WORKDIR /var/lib/minotour/apps

RUN git clone https://github.com/artic-network/artic-ncov2019.git

RUN git clone https://github.com/cov-lineages/pangolin.git

WORKDIR /var/lib/minotour/apps/artic-ncov2019

RUN conda env create -f environment.yml

WORKDIR /var/lib/minotour/apps/pangolin

RUN conda env create -f environment.yml

RUN conda install -y -c conda-forge uwsgi

RUN pip install -U setuptools

RUN pip install -U wheel

RUN pip install -r /var/lib/minotour/apps/minotourapp/requirements.txt && apt-get install -y nginx

COPY . /var/lib/minotour/apps/minotourapp/

WORKDIR /var/lib/minotour/apps/minotourapp

COPY docker/celeryd /etc/init.d/celeryd

COPY docker/celeryd.cnf /etc/default/celeryd

COPY docker/celerybeat /etc/init.d/celerybeat

COPY minotourapp_nginx.conf /etc/nginx/sites-available/

COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh

RUN ln -s /etc/nginx/sites-available/minotourapp_nginx.conf /etc/nginx/sites-enabled/

RUN cd extra/centrifuge-1.0.4-beta && make && cd ..

RUN useradd -ms /bin/bash celery

RUN usermod -a -G www-data root
RUN usermod -a -G root www-data

RUN python manage.py collectstatic

RUN chmod -R 755 extra/* && chmod 755 /etc/init.d/celeryd && chmod 640 /etc/default/celeryd && chmod 755 /etc/init.d/celerybeat
EXPOSE 8100

