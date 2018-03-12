#!/bin/bash

cd /opt/minotour/apps/minotourapp;

source ./envs.sh

/opt/minotour/envs/minotour/bin/python3 /opt/minotour/envs/minotour/bin/celery -A minotourapp worker -l info -B;
