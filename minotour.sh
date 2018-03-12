#!/bin/bash

cd /opt/minotour/apps/minotourapp;

source ./envs.sh

/opt/minotour/envs/minotour/bin/python3 /opt/minotour/envs/minotour/bin/gunicorn minotourapp.wsgi;
