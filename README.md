# Minotour

To do in this branch:
- Remove links from web `urls.py` and `views.py` that have been removed from the side menu
- Remove or repurpose `admin` dropdown menu
- Fix expanding table in active runs

---
minoTour webservice module


To run:

1. To serve the site (all run from root of site and each command must be run in own terminal):

python3 manage.py runserver port

2. For background tasks:

redis-server /usr/local/etc/redis.conf


celery -A minotourapp worker -l info -B

## Developer environment

This is a guideline for those creating development environments. Minotour
contains three main modules: 1) client, currently minFQ.py and minotour_api.py,
2) web, the interface that most end users have access and provides access to
active and archived runs; 3) rest api, that is a gateway connecting client and
web modules to data and core Minotour functionalities.

Minotour also makes use of a MySQL database, Celery (responsible for running
tasks), Redis (a database in memory similar to memcached, required by Celery),
Flower (optional - but a useful Celery task monitor tool).

In order to have a functional environment for development purposes, please
follow the instruction below.

* [Redis](https://redis.io/download) - Redis does not rely on the virtual
environment, so it can be installed anywhere.

* [MySQL Community edition](https://dev.mysql.com/downloads/) - Minotour
requires a MySQL server instance. It can run locally or in another
server.

* Python 3 - Minotour modules use Python 3, so make sure it is available.

* Virtual environment - Create a virtual environment to isolate Minotour
dependencies from other software installed on your computer.

```bash
mkdir ~/envs;

cd ~/envs;

python3 -m venv ~/envs/minotour;
```

* You do need to activate the virtual environment to run Minotour, Celery and
Flower, as well as, to install new dependencies.

```bash
source ~/envs/minotour/bin/activate;
```

* Clone Minotour repository from GitHub.

```bash
mkdir ~/projects;

cd ~/projects;

git clone https://github.com/py5gol/minotourws.git;

cd minotourws;

git checkout develop;
```

* Install Minotour dependencies.

```bash
pip install -r requirements.txt;
```

* Setup environment variables - many Minotour config parameters are defined in
environment, so that we do need to hardcode database users and password, or any
other information that is environment dependent. You can type it every time or
include in the ~/.bash_profile (please check if this is the correct file in
your environment).

```bash
export MT_DB_PASS=""
export MT_DB_PORT="3363"
export MT_DB_HOST="localhost"
export MT_DB_USER="root"
export MT_DB_NAME="minotourdb"
export MT_REFERENCE_LOCATION="~/data/minotour"
export MT_TWITCONSUMER_KEY=""
export MT_TWITTOKEN_SECRET=""
export MT_TWITCONSUMER_SECRET=""
export MT_TWITTOKEN=""
export MT_MAILGUN_ACCESS_KEY=""
export MT_MAILGUN_SERVER_NAME=""
```

* Now it is time to start the processes, so for each one, open a new terminal
and activate the virtual environment and follow the steps below.

** Start Celery

```bash

source ~/envs/minotour/bin/activate && cd ~/projects/minotourws/ && celery -A minotourapp worker -l info -B

```

** Start Flower

```bash
source ~/envs/minotour/bin/activate && cd ~/projects/minotourws/ && flow -A minotourapp --port=5555
```

** Start Minotour

```bash
source ~/envs/minotour/bin/activate && cd ~/projects/minotourws/ && pyth manage.py runserver 8100
```
