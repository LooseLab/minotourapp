FROM ubuntu:latest
MAINTAINER Rory Munro
SHELL [ "/bin/bash", "--login", "-c" ]
#!/bin/bash
ENV MT_DB_ENGINE='django.db.backends.mysql'
ENV MT_DB_HOST='localhost'
#ENV MT_DB_NAME='django_test'
ENV MT_DB_NAME='minotour'
ENV MT_DB_USER='dba'
ENV MT_DB_PASS=''
 # Mysql port. Defaults to 3306 at install.
ENV MT_DB_PORT='3306'
 # Boolean value. Run django in debug mode if True. NOT recommended, unless for development purposes.
ENV MT_DJANGO_DEBUG=0
# CURRENTLY UNUSED
ENV MT_MAILGUN_ACCESS_KEY=''
 # CURRENTLY UNUSED
ENV MT_MAILGUN_SERVER_NAME=''
ENV MT_MINIMAP2='/var/lib/minotour/apps/minotourapp/extra/minimap2-2.17_x64-linux/minimap2'
ENV MT_TWITCONSUMER_KEY='<Twitter consumer key>'
ENV MT_TWITCONSUMER_SECRET='<Twitter consumer secret>'
ENV MT_TWITTOKEN='<Twitter token>'
ENV MT_TWITTOKEN_SECRET='<Twitter token secret>'
# Generate a secret key with python manage.py shell -c 'from django.core.management import utils; print(utils.get_random_secret_key())'
ENV MT_SECRET_KEY='<Add secret key here>'
ENV MT_CENTRIFUGE="extra/centrifuge-1.0.4-beta/centrifuge"
 # Read centrifuge index from this directory
ENV MT_CENTRIFUGE_INDEX="/var/lib/minotour/apps/data"
 # Write Celery log file in this directory. Make sure you have write access
ENV MT_LOG_FOLDER='/var/lib/minotour/logs'
 # Redis url. Probably this value already, unless configured to run with sockets
ENV MT_CELERY_REDIS_URL='redis+socket:///redis_data/redis.sock'
 # Redis URl - where we store celery results
ENV MT_DJANGO_REDIS_URL='redis+socket:///redis_data/redis.sock'
# Number of Celery workers to start, no more than number of CPU cores
ENV MT_CELERY_WORKER_COUNT="2"
# Number of celery minimap2 threads to start, no more than num CPU cores
ENV MT_CELERY_MINIMAP2_WORKER_COUNT="1"
# Number of days since last flowcell activity until a flowcell is listed as inactive
ENV MT_TIME_UNTIL_INACTIVE="7"
# Number of days since last flowcell activity until a flowcell is archived and fastq data deleted
ENV MT_TIME_UNTIL_ARCHIVE="14"
# Boolean value of 1 or 0. If true (1) fastq data sent to the server is not stored in the database, and only aggregated metrics are kept
ENV MT_SKIP_SAVING_SEQUENCE="1"
# Boolean value, 1 or 0. If true, data can be uplaoded to the server. If false, it cannot be.
ENV MT_ALLOW_UPLOAD="1"
# Boolean value if True, failures in some celery tasks are tweeted to the administrator account for this minoTour instance.
ENV MT_TWEET_CELERY_FAILURES="1"
# if not setting up a server don't worry about this - it just shows the url on the client page for the example commands
ENV MT_SERVER_URL="localhost"
# The path and command used to activate conda in your system. On ubuntu this might be "/home/django/miniconda3/etc/profile.d/conda.sh" whereas on OSX it might be "/Users/username/.miniconda3/etc/profile.d/conda.sh"
ENV MT_CONDA_PREFIX="/root/miniconda3/etc/profile.d/conda.sh"
# by default these are stored in ../data/minimap2_indexes relative to the minotourapp base directory (not the app), but if specified indexes will be stored inside the given directory in a minimap2_indexes folder
ENV MT_MINIMAP2_INDEX_DIR="0"


##### Artic Task Environment Variables
## These variables are divided into those required for analysis and primer set configuration and those required for general configuration of the artic pipelines
# Read the Artic scheme directories, found where the Artic NCOV pipeline github repo was cloned. Example  /home/rory/Projects/artic-ncov2019/primer_schemes
ENV MT_ARTIC_SCHEME_DIR='/var/lib/minotour/data'
# The artic scheme name - e.g nCoV-2019 or SARS-CoV-2
ENV MT_ARTIC_SCHEME_NAME='SARS-CoV-2'
ENV MT_ARTIC_SCHEME_VER='V4'
 # The artic environment to activate to run artic commands.
ENV MT_ARTIC_ENV="artic"
 # If not setting up to use artic for Covid-19 this does not need to be set.
ENV MT_VoC_PATH="/var/lib/minotour/data"
ENV MT_ARTIC_MEDAKA_MODEL="r941_min_high_g351"
## These parameters specify the read lengths that minoTour will use to filter reads for the ARTIC pipeline. They are x
ENV MT_ARTIC_MIN_LEN=400
ENV MT_ARTIC_MAX_LEN=800
 # preferred coverage depth for artic
ENV MT_ARTIC_NORMALIZE=200
ENV MT_ARTIC_THREADS=4
ENV MT_ARTIC_MIN_AMPS_PERC="90"
ENV MT_IDEAL_READ_CONSTANT="100000"
# Boolean value if True, data that could potentially hold sensitive patient data that is fed into the artic pipeline is removed 12 hours after the run is finished.
ENV MT_DESTROY_ARTIC_EVIDENCE=1
ENV MT_COVERAGE_PER_AMPLICON="20"
 # Store the Artic results in this directory. Please ensure directory exists.
ENV MT_ARTIC_RESULTS_DIR='/var/lib/minotour/apps/data'
 # Maximum number of Artic tasks to be run at any one time
ENV MT_ARTIC_MAX_PIPELINES="1"
 # Number of hours until clearing artic files
ENV MT_ARTIC_TIME_UNTIL_CLEARING="12"
ENV TZ=Europe/London

ENV CONDA_DIR="/root/miniconda3"
RUN echo "$CONDA_DIR"
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV PYTHONUNBUFFERED 1

RUN mkdir -p /var/lib/minotour/apps/minotourapp /var/lib/minotour/logs /var/lib/minotour/data /var/lib/minotour/apps/data/references /var/lib/minotour/apps/minotour_static

WORKDIR /var/lib/minotour/apps/minotourapp

COPY requirements.txt /var/lib/minotour/apps/minotourapp/

RUN apt-get update && apt-get install -y git build-essential libz-dev python3-dev libssl-dev python3-pip mariadb-client default-libmysqlclient-dev

RUN apt-get install wget
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    chmod +x ~/miniconda.sh && \
    ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh
ENV PATH=$CONDA_DIR/bin:$PATH
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile
RUN conda init bash

WORKDIR /var/lib/minotour/apps

RUN git clone https://github.com/artic-network/fieldbioinformatics.git
RUN git clone https://github.com/cov-lineages/pangolin.git
WORKDIR /var/lib/minotour/apps/fieldbioinformatics
RUN conda install -c conda-forge mamba
RUN mamba env create -f environment.yml

WORKDIR /var/lib/minotour/apps/pangolin

RUN mamba env create -f environment.yml

RUN mamba install -y -c conda-forge uwsgi

RUN pip install -U setuptools

RUN pip install -U wheel

RUN pip install -r /var/lib/minotour/apps/minotourapp/requirements.txt && apt-get install -y nginx

COPY . /var/lib/minotour/apps/minotourapp/

WORKDIR /var/lib/minotour/apps/minotourapp
RUN mamba env create -f scripts/trees_environment.yml

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
RUN apt-get install -y build-essential python3-dev python3-pip python3-setuptools python3-wheel python3-cffi libcairo2 libpango-1.0-0 libpangocairo-1.0-0 libgdk-pixbuf2.0-0 libffi-dev shared-mime-info

RUN python manage.py collectstatic
RUN python scripts/init_ete3.py

RUN mkdir -p /redis_data

RUN chmod -R 755 extra/* && chmod 755 /etc/init.d/celeryd && chmod 640 /etc/default/celeryd && chmod 755 /etc/init.d/celerybeat && chmod 755 docker/init_conda.sh
EXPOSE 8100
ENTRYPOINT [ "/var/lib/minotour/apps/minotourapp/docker/init_conda.sh" ]

