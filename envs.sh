#!/bin/bash
export MT_DB_ENGINE='django.db.backends.mysql'
export MT_DB_HOST='localhost'
#export MT_DB_NAME='django_test'
export MT_DB_NAME='minotour' # Database name. Probably should be minotour
export MT_DB_USER='' # Database user. Could be root, recommend creating user for this, as per docs
export MT_DB_PASS='' # Database user password. Either root mysl password, or password of user that is created.
export MT_DB_PORT='3306' # Mysql port. Defaults to 3306 at install.
export MT_DJANGO_DEBUG=0 # Boolean value. Run django in debug mode if True. NOT recommended, unless for development purposes.
export MT_MAILGUN_ACCESS_KEY='' # CURRENTLY UNUSED
export MT_MAILGUN_SERVER_NAME='' # CURRENTLY UNUSED
export MT_MINIMAP2='</path/to/minotourapp/extra/minimap2-2.17_x64-linux/minimap2>' # use full path
export MT_REFERENCE_LOCATION='<path/to/references/directory/>' # Write uploaded references to this directory
export MT_TWITCONSUMER_KEY='<Twitter consumer key>'
export MT_TWITCONSUMER_SECRET='<Twitter consumer secret>'
export MT_TWITTOKEN='<Twitter token>'
export MT_TWITTOKEN_SECRET='<Twitter token secret>'
# Generate a secret key with python manage.py shell -c 'from django.core.management import utils; print(utils.get_random_secret_key())'
export MT_SECRET_KEY='<Add secret key here>'
export MT_CENTRIFUGE="</path/to/minotourapp/extra/centrifuge-1.0.4-beta/centrifuge>" # use full path
export MT_CENTRIFUGE_INDEX="/var/lib/minotour/apps/data" # Read centrifuge index from this directory
export MT_LOG_FOLDER='</path/to/logfile/directory/celery.log>' # Write Celery log file in this directory. Make sure you have write access
export MT_CELERY_BROKER_URL='redis://localhost:6379/0' # Redis url. Probably this value already, unless configured to run with sockets
export MT_CELERY_RESULT_BACKEND='redis://localhost:6379/0' # Redis URl - where we store celery results
export MT_CELERY_WORKER_COUNT="2" # Number of Celery workers to start, no more than number of CPU cores
export MT_CELERY_MINIMAP2_WORKER_COUNT="1" # Number of celery minimap2 threads to start, no more than num CPU cores
export MT_TIME_UNTIL_INACTIVE="7" # Number of days since last flowcell activity until a flowcell is listed as inactive
export MT_TIME_UNTIL_ARCHIVE="14" # Number of days since last flowcell activity until a flowcell is archived and fastq data deleted
export MT_SKIP_SAVING_SEQUENCE=0 # Boolean value of 1 or 0. If true (1) fastq data sent to the server is not stored in the database, and only aggregated metrics are kept
export MT_ALLOW_UPLOAD=1 # Boolean value, 1 or 0. If true, data can be uplaoded to the server. If false, it cannot be.
export MT_TWEET_CELERY_FAILURES=1 # Boolean value if True, failures in some celery tasks are tweeted to the administrator account for this minoTour instance.
export MT_SERVER_URL="localhost" # if not setting up a server don't worry about this - it just shows the url on the client page for the example commands

export MT_CONDA_PREFIX="" # The path and command used to activate conda in your system. On ubuntu this might be "/home/django/miniconda3/etc/profile.d/conda.sh" whereas on OSX it might be "/Users/username/.miniconda3/etc/profile.d/conda.sh"

export MT_MINIMAP2_INDEX_DIR="0" # by default these are stored in ../data/minimap2_indexes relative to the minotourapp base directory (not the app), but if specified indexes will be stored inside the given directory in a minimap2_indexes folder


##### Artic Task Environment Variables
## These variables are divided into those required for analysis and primer set configuration and those required for general configuration of the artic pipelines
export MT_ARTIC_SCHEME_DIR='<Path/To/Artic/Code/Scheme_Directories>' # Read the Artic scheme directories, found where the Artic NCOV pipeline github repo was cloned. Example  /home/rory/Projects/artic-ncov2019/primer_schemes
export MT_ARTIC_SCHEME_NAME='nCoV-2019' # The artic scheme name - e.g nCoV-2019 or
export MT_ARTIC_SCHEME_VER='V3' # The artic scheme version
export MT_ARTIC_ENV="artic-ncov2019" # The artic environment to activate to run artic commands.
export MT_VoC_PATH="<path/to/VoC/Definitions>" # If not setting up to use artic for Covid-19 this does not need to be set.
## These parameters specify the read lengths that minoTour will use to filter reads for the ARTIC pipeline. They are x
export MT_ARTIC_MIN_LEN=400
export MT_ARTIC_MAX_LEN=800
export MT_ARTIC_NORMALIZE=200 # preferred coverage depth for artic
export MT_ARTIC_THREADS=4 # Number of threads to run an artic job
export MT_ARTIC_MIN_AMPS_PERC="90"
export MT_IDEAL_READ_CONSTANT="100000"
export MT_DESTROY_ARTIC_EVIDENCE=0 # Boolean value if True, data that could potentially hold sensitive patient data that is fed into the artic pipeline is removed 12 hours after the run is finished.
export MT_COVERAGE_PER_AMPLICON="20"
export MT_ARTIC_RESULTS_DIR='/path/to/artic/directory' # Store the Artic results in this directory. Please ensure directory exists.
export MT_ARTIC_MAX_PIPELINES="1" # Maximum number of Artic tasks to be run at any one time
