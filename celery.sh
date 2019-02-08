#!/bin/bash

cd /var/lib/minotour/apps/minotourapp;

source ./envs.sh

/var/lib/minotour/envs/minotour/bin/python3 /var/lib/minotour/envs/minotour/bin/celery -A minotourapp worker -l info -B --uid 1016;

