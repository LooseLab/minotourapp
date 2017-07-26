# minotourws
minoTour webservice module


To run:

1. To serve the site (all run from root of site and each command must be run in own terminal):

python3 manage.py runserver port

2. For background tasks:

redis-server /usr/local/etc/redis.conf


celery -A minotourapp worker -l info -B


