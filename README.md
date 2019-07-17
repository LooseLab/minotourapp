# minoTour

## Warning - Work in progress

Minotour is a work in progress. Therefore features may change, and any bug reports/pull request would be gratefully recieved.

## Running minoTour with Docker

Firstly, make sure Docker is installed and running. You can check this running the command below:

```bash
docker --version
```

To following commands build and start several minoTour containers (database, Redis, celery, and web app).

```bash
git checkout https://github.com/LooseLab/minotourapp.git

cd minotourapp

docker-compose up
``` 

Now, you can access minoTour on http://localhost:10000.

## Development environment


To setup a local development environment, first clone or download and unzip the code from github found here - https://github.com/LooseLab/minotourapp.git

Checkout the development branch:

    git checkout develop

To run the python package mysqlclient, it is often necessary to have the following two dependencies installed, libmysqlclinet-dev and python3-dev. These can be installed with the following command:

    sudo apt-get install libmysqlclient-dev
    sudo apt-get install python3-dev

Create a virtual environment for the project dependencies and install them:

    cd /path/to/minotourapp/code
    python3 -m venv minotourenv
    source minotourenv/bin/activate
    pip install -r requirements.txt

Minotour contains three main modules:

1. The client, currently minFQ

2. web, the interface that most end users have access and provides access to active and archived runs

3. the rest api, that is a gateway connecting client and web modules to data and core Minotour functionalities.

Minotour also makes use of a MySQL database, Celery (responsible for running server tasks), Redis (a database in memory similar to memcached, required by Celery as a message broker), Flower (optional - but a useful Celery task monitor tool).


* [Redis](https://redis.io/download) - Minotour uses **Redis** as a cache system for the web module and also for Celery. Follow the instructions and make sure that the **redis-server** executable is available in PATH environment variable.

* [MySQL Community edition](https://dev.mysql.com/downloads/) - Minotour requires a MySQL server instance. It can run locally or on another server. Installing and configuring MySQL is not in the scope of this guide, but here is the official documentation https://dev.mysql.com/doc/mysql-getting-started/en/>`_ and another good tutorial can be found [here](https://www.digitalocean.com/community/tutorials/how-to-install-mysql-on-ubuntu-16-04). Once the server is up and running you can either choose to use the root user created during mySql initialisation in the below environmental variables, or create a user as follows in the mysql shell, logged in as the root user or an admin user:
    ```bash
    CREATE USER minotour;
    CREATE TABLE minotourdb;
    GRANT ALL PRIVILEGES ON minotourdb TO 'minotour'@'localhost' IDENTIFIED BY '<password>'

    ```
    

* [Minimap2](https://github.com/lh3/minimap2) - Minotour uses **Minimap2** to run fast alignment, and to do metagenomics target validation.

* [Python 3](https://www.python.org) - Minotour uses **Python >=3.6**, so make sure it is available on your system.

* You will need to have the virtual environment activated to run Minotour, Celery and Flower, as well as to install the new dependencies.

* Create Minotour data directory - currently, Minotour uses this folder to keep track of the genome references available:

    mkdir -p ~/data/minotour;

---------------------
Metagenomics analyses
---------------------

minoTour uses [Centrifuge](https://ccb.jhu.edu/software/centrifuge/) to run metagenomics analyses.
Before creating new metagenomics tasks, there are a few requirements that need to be completed.

- Centrifuge index - Choose one of the centrifuge indexes (ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data) and save in the local disk. We recommend the compressed index p_compressed

- Centrifuge application - [Download](https://github.com/infphilo/centrifuge/releases) and compile the most recent version.

- Set the environment variables MT_CENTRIFUGE (the path to the executable) and MT_CENTRIFUGE_INDEX (the path to indexes without the suffix) as below:

    export MT_CENTRIFUGE="/home/user/centrifuge-1.0.4-beta/centrifuge"

    export MT_CENTRIFUGE_INDEX="/home/user/centrifuge_indexes/p_compressed"

    The MT_CENTRIFUGE_INDEX above points to the following indexes:
    
    p_compressed.1.cf
    
    p_compressed.2.cf
    
    p_compressed.3.cf
    
    p_compressed.4.cf

- minoTour uses the ete3 package, that needs access to the internet to download the NCBI species database.

To force the download (this step just need to be executed once), type the following instructions into the python interpreter:
```bash
    import ete3

    from ete3 import NCBITaxa

    n = NCBITaxa()
```

If you are operating on a server, you may wish to raise the priority of the centrifuge process, as otherwise it may not be allocated free RAM in a workable timeframe.

In order to do this you'll need to edit the /etc/security/limits.conf file to allow non root users to set a processes niceness index to less than 0.

To do this you will need root permissions. Run the following terminal commands:

```bash
sudo vim /etc/security/limits.conf
``` 
Which will show you a that looks like this:
```bash
    #
#/etc/security/limits.conf
#
#Each line describes a limit for a user in the form:
#
#<domain>        <type>  <item>  <value>
#
#Where:
#<domain> can be:
#        - a user name
#        - a group name, with @group syntax
#        - the wildcard *, for default entry
#        - the wildcard %, can be also used with %group syntax,
#                 for maxlogin limit
#        - NOTE: group and wildcard limits are not applied to root.
#          To apply a limit to the root user, <domain> must be
#          the literal username root.
#
#<type> can have the two values:
#        - "soft" for enforcing the soft limits
#        - "hard" for enforcing hard limits
#
#<item> can be one of the following:
#        - core - limits the core file size (KB)
#        - data - max data size (KB)
#        - fsize - maximum filesize (KB)
#        - memlock - max locked-in-memory address space (KB)
#        - nofile - max number of open files
#        - rss - max resident set size (KB)
#        - stack - max stack size (KB)
#        - cpu - max CPU time (MIN)
#        - nproc - max number of processes
#        - as - address space limit (KB)
#        - maxlogins - max number of logins for this user
#        - maxsyslogins - max number of logins on the system
#        - priority - the priority to run user process with
#        - locks - max number of file locks the user can hold
#        - sigpending - max number of pending signals
#        - msgqueue - max memory used by POSIX message queues (bytes)
#        - nice - max nice priority allowed to raise to values: [-20, 19]
#        - rtprio - max realtime priority
#        - chroot - change root to directory (Debian-specific)
#
#<domain>      <type>  <item>         <value>
#

#*               soft    core            0
#root            hard    core            100000
#*               hard    rss             10000
#@student        hard    nproc           20
#@faculty        soft    nproc           20
#@faculty        hard    nproc           50
#ftp             hard    nproc           0
```

Add the following line to the file, which will allow non root users to set a limit of -10 to the niceness index, giving the process a higher priority.

```bash
    *	-	nice	-10
```

This step isn't necessary, but it might be useful if the centrifuge index is taking a long time to load. 


    

--------------------------------
Environmental config and running
--------------------------------

* Setup environment variables - many Minotour config parameters are stored in the environment, so that we do need to hardcode database users and password, or any other information that is environment dependent. You can create a bash script file, for example envs.sh, place it in the main application directory (The directory that contains manage.py) and copy the following into it, or include them in the ~/.bash_profile or ~/.bash_rc (please check if this is the correct file in your environment).::
    ```bash
        #!/bin/bash
        export MT_DB_ENGINE='django.db.backends.mysql'
        export MT_DB_HOST='localhost'
        export MT_DB_NAME='minotourdb'
        export MT_DB_PASS='<minotourdb password>'
        export MT_DB_PORT='3306'
        export MT_DB_USER='minotour'
        export MT_DJANGO_DEBUG='True'
        export MT_MAILGUN_ACCESS_KEY=''
        export MT_MAILGUN_SERVER_NAME=''
        export MT_MINIMAP2='/usr/bin/minimap2'
        export MT_REFERENCE_LOCATION='/home/rory/data/minotour'
        export MT_SECRET_KEY=''
        export MT_TWITCONSUMER_KEY=''
        export MT_TWITCONSUMER_SECRET=''
        export MT_CENTRIFUGE="</path/to/centrifuge/executable>"
        export MT_CENTRIFUGE_INDEX="/path/to/centrifuge/indexes"
        export MT_TWITTOKEN=''
        export MT_TWITTOKEN_SECRET=''
        export MT_LOG_FOLDER='/path/to/where/you/want/logs'
        export MT_CELERY_BROKER_URL='redis://localhost:6379/0'
        export MT_CELERY_RESULT_BACKEND='redis://localhost:6379/0'
    ```
    
* Now it is time to start the processes, we suggest opening a new terminal for each command. If you chose to create the environmental variable bash file, add the following to the beginning of the celery, flower and Minotour commands to set the environmental variables::
    ```bash
      . envs.sh &&
    ```

* MySQL - make sure it is running **AND the database was created** (check the docs mentioned above).

* Redis:
    ```bash
      redis-server &
    ```

* Create tables and administrator account:
    ```bash
        cd /path/to/minotour/code;
        source ~/minotourenv/bin/activate;
        python3 manage.py makemigrations;
        python3 manage.py migrate;
        python3 manage.py loaddata fixtures/auxiliary_data.json
        python3 manage.py createsuperuser
    ```

* Start Celery::
    ```bash
        cd /path/to/minotour/code/ && source minotourenv/bin/activate && celery -A minotourapp worker -l info -B
    ```

* Start Flower::
    ```bash
        cd /path/to/minotour/code/ && source minotourenv/bin/activate && && flower -A minotourapp --port=5555
    ```

* Start Minotour::
    ```bash
        cd /path/to/minotour/code
        source minotourenv/bin/activate
        python manage.py runserver 8100
    ```

* Time to test - if everything worked well, you should be able to access the web interface on http://localhost:8100.








