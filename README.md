# Minotour

To do in this branch:
- Remove links from web `urls.py` and `views.py` that have been removed from the side menu
- Remove or repurpose `admin` dropdown menu
- Fix expanding table in active runs

---
minoTour webservice module


To run (we hope):

* To serve the site (all run from root of site and each command must be run in own terminal):

```bash
python3 manage.py runserver port
```

* For background tasks:

```bash
redis-server /usr/local/etc/redis.conf

celery -A minotourapp worker -l info -B

flow -A minotourapp --port=5555
```

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

* [Redis](https://redis.io/download) - Minotour uses **Redis** as a cache system
for the web module and also for Celery. Follow the instructions and make sure
that the **redis-server** executable is available in PATH environment variable.

* [MySQL Community edition](https://dev.mysql.com/downloads/) - Minotour
requires a MySQL server instance. It can run locally or in another
server. Installing and configuring MySQL is not in the scope of this guideline, but here is official
documentation [Link](https://dev.mysql.com/doc/mysql-getting-started/en/) and
another good tutorial is [Link](https://www.digitalocean.com/community/tutorials/how-to-install-mysql-on-ubuntu-16-04).

* [Minimap2](https://github.com/lh3/minimap2) - Minotour uses **Minimap2** to run
fast alignment, so you are not running alignments, just skip this step.
Otherwise, follow the instructions and make sure that the **minimap2** executable
is available in PATH environment variable.

* [Bwa](https://github.com/lh3/bwa) - Minotour uses **bwa** to run
alignment, so you are not running alignments, just skip this step.
Otherwise, follow the instructions and make sure that the **bwa** executable
is available in PATH environment variable.

* [Python 3](https://www.python.org) - Minotour uses Python 3, so make sure it
is available in your system.

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

* Create Minotour data directory - currently, Minotour uses this folder to
keep track of the genome references available.

```bash
mkdir -p ~/data/minotour;
```

* Setup environment variables - many Minotour config parameters are defined in
environment, so that we do need to hardcode database users and password, or any
other information that is environment dependent. You can type it every time or
include in the ~/.bash_profile or ~/.bash_rc (please check if this is the correct file in
your environment).

```bash
export MT_DB_PASS="<your must inform the root or other db user password>"
export MT_DB_PORT="3363"
export MT_DB_HOST="localhost"
export MT_DB_USER="root or other db user"
export MT_DB_NAME="minotourdb"
export MT_REFERENCE_LOCATION="~/data/minotour"
export MT_TWITCONSUMER_KEY=""
export MT_TWITTOKEN_SECRET=""
export MT_TWITCONSUMER_SECRET=""
export MT_TWITTOKEN=""
export MT_MAILGUN_ACCESS_KEY=""
export MT_MAILGUN_SERVER_NAME=""
```

* Now it is time to start the processes, and we suggest opening a new terminal
for each command.

    * MySQL - make sure it is running **AND the database was created** (check the docs mentioned above).

    * Redis

    ```bash
    redis-server &
    ```
    * Create tables

    ```bash
    source ~/envs/minotour/bin/activate;

    cd ~/projects/minotourws;

    python manage.py migrate;
    ```
    
    * Start Celery

    ```bash
    source ~/envs/minotour/bin/activate && cd ~/projects/minotourws/ && celery -A minotourapp worker -l info -B
    ```

    * Start Flower

    ```bash
    source ~/envs/minotour/bin/activate && cd ~/projects/minotourws/ && flower -A minotourapp --port=5555
    ```

    * Start Minotour

    ```bash
    source ~/envs/minotour/bin/activate;

    cd ~/projects/minotourws;

    python manage.py runserver 8100;
    ```

* Time to test - if everything worked well, you should be able to access the web
interface on http://localhost:8100. Login with user **demo** and password **asdf1234**.
