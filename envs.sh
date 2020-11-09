#!/bin/bash
export MT_DB_ENGINE='django.db.backends.mysql'
export MT_DB_HOST='localhost'
#export MT_DB_NAME='django_test'
export MT_DB_NAME='minotour'
export MT_DB_USER='dba>'
export MT_DB_PASS='test'
export MT_DB_PORT='3306'
export MT_DJANGO_DEBUG=0
export MT_MAILGUN_ACCESS_KEY=''
export MT_MAILGUN_SERVER_NAME=''
export MT_MINIMAP2='extra/minimap2-2.17_x64-linux/minimap2'
export MT_REFERENCE_LOCATION='<path/to/references/directory/>' # Write uploaded references to this directory
export MT_TWITCONSUMER_KEY='<Twitter consumer key>'
export MT_TWITCONSUMER_SECRET='<Twitter consumer secret>'
export MT_TWITTOKEN='<Twitter token>'
export MT_TWITTOKEN_SECRET='<Twitter token secret>'
# Generate a secret key with python manage.py shell -c 'from django.core.management import utils; print(utils.get_random_secret_key())'
export MT_SECRET_KEY='<Add secret key here>'
export MT_CENTRIFUGE="extra/centrifuge-1.0.4-beta/centrifuge"
export MT_CENTRIFUGE_INDEX="/var/lib/minotour/apps/data" # Read centrifuge index from this directory
export MT_LOG_FOLDER='/var/lib/minotour/logs' # Write Celery log file in this directory
export MT_CELERY_BROKER_URL='redis://localhost:6379/0'
export MT_CELERY_RESULT_BACKEND='redis://localhost:6379/0'
#export MT_ARTIC_RESULTS_DIR='/media/rory/Bioinformatics'
export MT_ARTIC_RESULTS_DIR='/var/lib/minotour/apps/data' # Store the Artic results in this directory
export MT_ARTIC_SCEHEME_DIR='<Path/To/Artic/Code/Scheme_Directories>' # Read the Artic scheme directories, found where the Artic NCOV piprline github repo was clones. Example  Projects/artic-ncov2019/primer_schemes
export MT_ARTIC_MAX_PIPELINES="1" # Maximum number of Artic tasks to be run at any one time
export MT_CELERY_WORKER_COUNT="2" # Number of Celery workers to start, no more than number of CPU cores
export MT_CELERY_MINIMAP2_WORKER_COUNT="1" # Number of celery minimap2 threads to start, no more than num CPU cores
export MT_TIME_UNTIL_INACTIVE="7" # Number of days since last flowcell activity until a flowcell is listed as inactive
export MT_TIME_UNTIL_ARCHIVE="14" # Number of days since last flowcell activity until a flowcell is archived and fastq data deleted
export MT_SKIP_SAVING_SEQUENCE=0
