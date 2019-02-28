#!/bin/sh

# pip3 install pip --upgrade;
pip3 install -r requirements.txt;

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
export MT_MINIMAP2="/var/lib/minotour/apps/minotourapp/extra/minimap2/minimap2";
export MT_CELERY_BROKER_URL='redis://redis-minotour-instance:6379/0'
export MT_CELERY_RESULT_BACKEND='redis://redis-minotour-instance:6379/0'
export MT_CENTRIFUGE="/var/lib/minotour/apps/minotourapp/extra/centrifuge-1.0.4-beta/centrifuge"
export MT_CENTRIFUGE_INDEX="/var/lib/minotour/data/p_compressed"

celery -A minotourapp worker -l info;
