#!/bin/bash

#/etc/init.d/celeryd start
#/etc/init.d/celerybeat start
uwsgi --ini minotourapp_uwgsi.ini
nginx -g 'daemon off;'
