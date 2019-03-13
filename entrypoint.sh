#!/bin/sh
pip install pip --upgrade;
pip install -r requirements.txt;
export MT_DB_ENGINE="django.db.backends.mysql";
export MT_DB_NAME="minotour";
export MT_DB_USER="root";
export MT_DB_PASS="test";
export MT_DB_PORT="3306";
export MT_DB_HOST="db-minotour-instance";
export MT_SECRET_KEY='';
export MT_DJANGO_DEBUG='True';
export MT_MAILGUN_ACCESS_KEY='';
export MT_MAILGUN_SERVER_NAME='';
export MT_TWITTOKEN='';
export MT_TWITTOKEN_SECRET='';
export MT_TWITCONSUMER_KEY='';
export MT_TWITCONSUMER_SECRET='';
export MT_REFERENCE_LOCATION='/var/lib/minotour/data';
export MT_MINIMAP2="/var/lib/minotour/minimap2/minimap2";
export MT_CELERY_BROKER_URL='redis://redis-minotour-instance:6379/0'
export MT_CELERY_RESULT_BACKEND='redis://redis-minotour-instance:6379/0'
python3 manage.py makemigrations;
python3 manage.py migrate;
python3 manage.py loaddata fixtures/auxiliary_data.json;
exec python3 manage.py runserver 0.0.0.0:8000;
