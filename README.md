# minoTour

Minotour is a work in progress, therefore features may change, and any bug reports/feature requests and pull requests will be gratefully recieved.

## Running minoTour with Docker
minoTour is developed on docker so docker must be installed and rightly set up to develop and test miniTour

Docker is a set of platform as a service products that use OS-level virtualization to deliver software in packages called containers. 
Containers are isolated from one another and bundle their own software, libraries and configuration files;
they can communicate with each other through well-defined channels.

If Docker is not installed - instructions for installing can be found here <https://docs.docker.com/install/>
Make sure Docker is installed and running at the end of the installation. You can check this by running the command below:

```bash
docker --version
```
## Building and starting Docker containers
The following commands build and start several minoTour docker containers and create an admin user.
clone the minoTour repository unto your local machine 
```bash
git clone https://github.com/LooseLab/minotourapp.git
``` 

Change to the minotourapp directory and then use the git checkout command to switch to the develop branch
```bash
cd minotourapp

git checkout develop
``` 
 
Use the following commands to build and start docker services 
- <!-- Use 'docker-compose up web' if celere is not running -->-

```bash
docker-compose build

docker-compose up -d 
``` 

Create an admin user:

```bash
docker-compose exec web python manage.py createsuperuser
``` 

After a few seconds, you should be able to access minoTour at <http://localhost:8100/>


### Uploading references and other files

Once a instance is up and running, a local data folder should have appeared alongside the minotourapp directory. This is for any reference files and the centrifuge index (if you want to run real).
 We do not distribute this index, it can be dowloaded from [JHU](ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed_2018_4_15.tar.gz).
 
 After placing any files into the localdata folder on you host system, run the following the commands to
 make them available to minoTour.
 
 ```bash 
docker exec celery_worker python3 manage.py add_reference /var/lib/minotour/apps/localdata
```

To add the centrifuge index to the docker container, create a folder called centrifuge_files in the localdata folder on your host system.
Place the 4 .cf files in the newly created folder. 
 

## Development environment


## Setting up development environment

To set up a development environment, first clone or download the code from [our github page](https://github.com/LooseLab/minotourapp.git "Looselab's github page").

Checkout the development branch:

    git checkout develop

To use the python package mysqlclient, it is necessary to have two dependencies installed, libmysqlclinet-dev and python3-dev. These can be installed with the following command:

    sudo apt-get install libmysqlclient-dev
    sudo apt-get install python3-dev

Now create a virtual environment for the project dependencies and install them:

    cd /path/to/minotourapp/code
    python3 -m venv minotourenv
    source minotourenv/bin/activate
    pip install -r requirements.txt

minoTour is seperated into three main modules:

1. The upload client, minFQ

2. The minoTour web app, the interface that provides the end user access to active and archived runs

3. The rest API, that connects minFQ and the web app to the back end data and analyses.

Minotour uses a MySQL database, Celery (asynchronous workers responsible for running server tasks), Redis (a database in memory similar to memcached, required by Celery as a message broker), Flower (optional - a useful Celery task monitor tool).


* [Redis](https://redis.io/download) - Minotour uses **Redis** as a cache memory system and message broker for Celery. Follow the linked installation instructions.

* [MySQL Community edition](https://dev.mysql.com/downloads/) - Minotour requires a running instance of mySQL server. To set up mySQL community edition follow the instructions found in the [official installation documentation](https://dev.mysql.com/doc/mysql-getting-started/en/). Another good tutorial can be found [here](https://www.digitalocean.com/community/tutorials/how-to-install-mysql-on-ubuntu-16-04). Once the server is up and running, you can either use the root user created during set up, setting it in the below environment variables file.   
 Otherwise you can create a user using the mysql shell, whilst logged in as the root user or as an admin user:
    ```bash
    CREATE USER minotour;
    CREATE TABLE minotourdb;
    GRANT ALL PRIVILEGES ON minotourdb TO 'minotour'@'localhost' IDENTIFIED BY '<password>'

    ```
    

* [Minimap2](https://github.com/lh3/minimap2) - Minotour uses **Minimap2** to run alignments, and for metagenomics target validation.

* [Python 3](https://www.python.org) - Minotour uses **Python >=3.6**, so make sure it is available on your system.



## Conda Environments for Artic analyses
To run Artic analyses, conda must be installed, and the Artic and Pangolin environments must be installed in it. We recommend installing miniconda to avoid bloat.

The miniconda installation instructions can be found [here](https://docs.conda.io/en/latest/miniconda.html)

The [Artic](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html) instructions, and [Pangolin](https://github.com/cov-lineages/pangolin#install-pangolin) instructions.


---------------------
Metagenomics analyses
---------------------

minoTour uses [Centrifuge](https://ccb.jhu.edu/software/centrifuge/) to run metagenomics analyses.
Before creating a metagenomics tasks, there are a few requirements that need to be met.

- A centrifuge index is required - Choose one of the centrifuge [indexes](ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data), download it and save it to the local disk. We recommend the compressed index, p_compressed. This index contains Bacteria and Archaea down to the species level.

- Centrifuge itself - [Download](https://github.com/infphilo/centrifuge/releases) and compile the most recent version.

- Set the environment variables MT_CENTRIFUGE (the path to the executable) and MT_CENTRIFUGE_INDEX (the path to indexes without the suffix) to the location of your chosen Index and the centrifuge executable:

    export MT_CENTRIFUGE="/path/to/centrifuge-executable/centrifuge-1.0.4-beta/centrifuge"

    export MT_CENTRIFUGE_INDEX="/path/to/centrifuge_indexes/p_compressed"

    The MT_CENTRIFUGE_INDEX above points to the following indexes:
    
    * p_compressed.1.cf
    
    * p_compressed.2.cf
    
    * p_compressed.3.cf
    
    * p_compressed.4.cf

- minoTour uses the ete3 package for Id'ing species from there taxonomic IDs. The installation of ete3's database during the first Metagenomics analysis will happen automatically, but requires internet access.

 - Alternatively, to force the download (this step just need to be executed once), type the following instructions into the python shell in the minoTour virtual environment:
```bash
    import ete3

    from ete3 import NCBITaxa

    n = NCBITaxa()
```

If you are operating on a server, you may wish to raise the priority of the centrifuge process, as otherwise it may not be allocated enough RAM to load the whole index into memory.

This step isn't necessary, but it might be useful if the centrifuge index is taking a long time to load. 

In order to do this you'll need to edit the /etc/security/limits.conf file to allow non root users to set a processes niceness index to less than 0.

To do this you will need root permissions. Run the following terminal commands:

```bash
sudo vim /etc/security/limits.conf
``` 
Which will show you a file that looks like this:
```bash
    
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

Add the following line to the bottom of the file, which will allow non root users to set a limit of -10 to the niceness index, giving the process a higher priority.

```bash
    *	-	nice	-10
```



    

--------------------------------
Environment config
--------------------------------

* Setting up environment variables - many of minoTour's config parameters are stored in the environment. To store these variables,
 create a bash script file, for example envs.sh, place it in the main application directory (The directory that contains manage.py)
  and copy the following into it:
    ```bash
        #!/bin/bash
        export MT_DB_ENGINE='django.db.backends.mysql'
        export MT_DB_HOST='localhost'
        export MT_DB_NAME='<Database name>'
        export MT_DB_PASS='<minotourdb password>'
        export MT_DB_PORT='3306'
        export MT_DB_USER='minotour'
        export MT_DJANGO_DEBUG='True'
        export MT_MAILGUN_ACCESS_KEY=''
        export MT_MAILGUN_SERVER_NAME=''
        export MT_MINIMAP2='</Path/To/Minimap2/Executable>'
        export MT_REFERENCE_LOCATION='<Path/To/Reference/Directory>'
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
    
-------------------
Starting the server
-------------------
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
    
* If you created the environment variable bash file, add the following to the beginning of the celery, flower and minoTour commands to set the environment variables:
    ```bash
      . envs.sh &&
    ```

* Start Celery::
    ```bash
        cd /path/to/minotour/code/ 
        source minotourenv/bin/activate
        celery -A minotourapp worker -l info -B
    ```

* Start Flower::
    ```bash
        cd /path/to/minotour/code/
        source minotourenv/bin/activate
        flower -A minotourapp --port=5555
    ```

* Start Minotour::
    ```bash
        cd /path/to/minotour/code
        source minotourenv/bin/activate
        python manage.py runserver 8100
    ```

* Time to test - if all went well, you should be able to access the web interface on http://localhost:8100.
