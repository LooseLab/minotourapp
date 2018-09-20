# minoTour



### minoTour webservice module


To run:

* To serve the site in production (run from root directory of Django Project, all background tasks must be running):

```bash
python3 manage.py runserver port
```

* For background tasks:

```bash
redis-server /usr/local/etc/redis.conf

celery -A minotourapp worker -l info -B

flower -A minotourapp --port=5555
```

## Developer environment

This is a guideline for those creating development environments. Minotour
contains three main modules: 1) client, currently minFQ.py and minotour_api.py,
2) web, the interface that most end users have access and provides access to
active and archived runs; 3) the rest api, that is a gateway connecting client and
web modules to data and core Minotour functionalities.

Minotour also makes use of a MySQL database, Celery (responsible for running
server tasks), Redis (a database in memory similar to memcached, required by Celery as a message broker),
Flower (optional - but a useful Celery task monitor tool).

In order to have a functional environment for development purposes, please
follow the instruction below.

* [Redis](https://redis.io/download) - Minotour uses **Redis** as a cache system
for the web module and also for Celery. Follow the instructions and make sure
that the **redis-server** executable is available in PATH environment variable.

* [MySQL Community edition](https://dev.mysql.com/downloads/) - Minotour
requires a MySQL server instance. It can run locally or on another
server. Installing and configuring MySQL is not in the scope of this guide, but here is the official
documentation [Link](https://dev.mysql.com/doc/mysql-getting-started/en/) and
another good tutorial can be found here [Link](https://www.digitalocean.com/community/tutorials/how-to-install-mysql-on-ubuntu-16-04).
Once the server is up and running you can either choose to use the root user created during mySql initialisation in the below environmental variables, or create a user as follows in the mysql shell, logged in as the root user or an admin user:

```bash
CREATE USER minotour;
CREATE TABLE minotourdb;
GRANT ALL PRIVILEGES ON minotourdb TO 'minotour'@'localhost' IDENTIFIED BY '<password>'
```

* [Minimap2](https://github.com/lh3/minimap2) - Minotour uses **Minimap2** to run
fast alignment, so if you are not running alignments, just skip this step.
Otherwise, follow the instructions and make sure that the **minimap2** executable
is available in the PATH environment variable.

* [Python 3](https://www.python.org) - Minotour uses Python 3, so make sure it
is available on your system.

* Virtual environment - Create a virtual environment to isolate Minotour
dependencies from other software installed on your computer.

```bash
mkdir ~/envs;

cd ~/envs;

python3 -m venv ~/envs/minotour;
```

* You will need to activate the virtual environment to run Minotour, Celery and
Flower, as well as to install the new dependencies.

```bash
source ~/envs/minotour/bin/activate;
```

* Clone Minotour repository from GitHub.

```bash
mkdir ~/projects;

cd ~/projects;

git clone https://github.com/LooseLab/minotourapp.git

cd minotourws;

git checkout develop;
```

* To run the python package mysqlclient, it is often necessary to have the following
two dependencies installed, libmysqlclinet-dev and python3-dev. These can be installed with the following command
```bash
sudo apt-get install libmysqlclient-dev

sudo apt-get install python3-dev
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
other information that is environment dependent. You can create a bash script file, for example envs.sh, place it in the main application directory (The directory that contains manage.py) and
copy the following into it, or include in the ~/.bash_profile or ~/.bash_rc (please check if this is the correct file in
your environment).


```bash
export MT_DB_PASS="<mySQL users password>"
export MT_DB_PORT="3363"
export MT_DB_HOST="localhost"
export MT_DB_USER="root or other db user, i.e minotour created above"
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
for each command. If you chose to create the environmental variable bash file, add the following to the beginning of the celery, flower and Minotour commands to set the environmental variables - 

```bash
. envs.sh &&
```

    * MySQL - make sure it is running **AND the database was created** (check the docs mentioned above).

    * Redis

    ```bash
    redis-server &
    ```
    * Create tables

    ```bash
    source ~/envs/minotour/bin/activate;

    cd ~/projects/minotourws;
    python manage.py makemigrations;
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
interface on http://localhost:8100.

## Centrifuge

"Specify the "
