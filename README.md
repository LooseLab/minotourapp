# minoTour

minoTour is a **work in progress**, therefore features may change, and any bug reports/feature requests and pull requests will be gratefully recieved.

minoTour is a web-based real-time laboratory information management system (LIMS) for Oxford Nanopore Technology (ONT) sequencers built using the Django framework. minoTour can be hosted offline and locally on a laptop, or online and accessed remotely.  It aims to enable monitoring of Nanopore sequencers in real time, allowing for analysis of performance and basic downstream data analysis. minoTour can monitor the activity of a sequencer in real time independent of analysing basecalled files providing a breakdown and analysis of live sequencing metrics via integration with ONT’s minKNOW API and parsing of sequence files as they are generated. These features enable remote control of a sequencer, as well as tracking flow cell QC and all messages received by MinKNOW during a sequencing run. The client for monitoring and uploading data, minFQ, is available via PyPI. Simple analyses, such as base calling summaries, alignment and metagenomics can be completed in real-time or after sequencing. minoTour is open source (https://www.github.com/LooseLab/minotourapp)

# Table of contents
* [Preconfiguration](#preconf)  
* [Setting up minoTour](#setup)  
    * [Set up minoTour and requirements](#set-up-minotour)
    * [Set up environment variables](#setting-up-the-environment-variable-file)
    * [Create database tables and admin](#create-database-tables-and-administrator-account)
    * [Conda environments](#conda-environments-for-artic-analyses)
* [Running minotour with tmux](#running-minotour-with-tmux-recommended)
* [Starting the server manually](#starting-the-server-manually)
* [Final configurations](#final-configurations)
    * [Metagenomics analyses](#metagenomics-analyses)
    * [Raise allowed Niceness index](#raise-the-niceness-index-allowed)
    * [Setting up a twitter API account](#setting-up-twitter)

## <a name="preconf"></a>Preconfiguration

Minotour uses MySQL as a database backend, Celery (asynchronous workers responsible for running server tasks), Redis (a database in memory similar to memcached, required by Celery as a message broker), Flower (optional - a useful Celery task monitor tool).

The following tools must be installed:

* [Redis](https://redis.io/download) - Minotour uses **Redis** as a memory cache system and message broker for Celery. Follow the linked installation instructions.

* [MySQL Community edition](https://dev.mysql.com/downloads/) - Minotour requires a running instance of mySQL server. To set up mySQL community edition follow the instructions found in the [official installation documentation](https://dev.mysql.com/doc/mysql-getting-started/en/). Another good tutorial can be found [here](https://www.digitalocean.com/community/tutorials/how-to-install-mysql-on-ubuntu-16-04). Once the server is up and running, you can either use the root user created during set up, setting it in the below environment variables file.   
 Otherwise you can create a user using the mysql shell, whilst logged in as the root user or as an admin user:
    ```bash
    CREATE USER minotour;
    CREATE TABLE minotourdb;
    GRANT ALL PRIVILEGES ON minotourdb TO 'minotour'@'localhost' IDENTIFIED BY '<password>'
    ```

* [Python 3](https://www.python.org) - Minotour uses **Python >=3.6**, so please make sure it is available on your system.

To upload data, it is also necessary to install the upload client, [minFQ](https://github.com/LooseLab/minotourcli)
##<a name="setup"></a> Setting up the minoTour environment


### Set up minoTour
To set up a development environment, first clone or download the code from [our github page](https://github.com/LooseLab/minotourapp.git "Looselab's github page").

Checkout the development branch:

    git checkout develop

To use the python package mysqlclient, it is necessary to have two dependencies installed, libmysqlclinet-dev and python3-dev. These can be installed with the following command:

    sudo apt-get install libmysqlclient-dev
    sudo apt-get install python3-dev

Now create a virtual environment for the project dependencies and install them:

    cd /path/to/minotourapp/
    python3 -m venv minotourenv
    source minotourenv/bin/activate
    pip install -r requirements.txt
    
### Setting up the environment variable file

Setting up environment variables - many of minoTour's config parameters are stored in the environment. There is a file to be configured,
[envs.sh](envs.sh)
For the secret key, we recommend generating your own, which can be done as follows:
```bash
$ python manage.py shell -c 'from django.core.management import utils; print(utils.get_random_secret_key())'
```
And stored under MT_SECRET_KEY in envs.sh.

### Create database tables and administrator account:
Once the environment variable file has been configured, we need to make the tables to store the data and if desired we need to create and admin account to manage the data. Luckily django makes this easy...
```bash
    cd /path/to/minotour/code;
    source ~/minotourenv/bin/activate;
    source envs.sh
    python3 manage.py makemigrations;
    python3 manage.py migrate;
    python3 manage.py loaddata fixtures/auxiliary_data.json
    python3 manage.py createsuperuser
```

## Conda Environments for Artic analyses
To run Artic analyses, conda must be installed, and the Artic and Pangolin environments must be installed in it. We recommend installing miniconda to avoid bloat.

The miniconda installation instructions can be found [here](https://docs.conda.io/en/latest/miniconda.html)

The [Artic](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html) instructions, and [Pangolin](https://github.com/cov-lineages/pangolin#install-pangolin) instructions.

## Running minoTour with tmux (Recommended)

tmux is a [terminal multiplexer](https://en.wikipedia.org/wiki/Tmux) for Unix-like operating systems, namely linux and mac OSX. Installation instructions can be found [here for linux](https://tmuxguide.readthedocs.io/en/latest/tmux/tmux.html).
To install tmux on mac, use homebrew:
```bash
brew install tmux
```
Once installed, minoTour has a [bash script](scripts/run_minotour_tmux.sh) that runs all the tmux commands to open minoTour in Tmux. The configuration of the file is necessary to provide all the correct paths.
It can be run simply from the **main minotourapp directory**:
```bash
./scripts/run_minotour_tmux.sh
```

## Starting the server manually

We recommend using the tmux script for starting the server, but if you wish you can start each if the below processes in terminal.
* Redis:
    ```bash
      redis-server &
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


## Final configurations

---------------------
Metagenomics analyses
---------------------

minoTour uses [Centrifuge](https://ccb.jhu.edu/software/centrifuge/) to run metagenomics analyses.
Before creating a metagenomics tasks, there are a few requirements that need to be met.

- A centrifuge index is required - Choose one of the centrifuge [indexes](ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data), download it and save it to the local disk. We recommend the compressed index, p_compressed. This index contains Bacteria and Archaea down to the species level.
Instructions to build a custom index can be found [here](http://www.ccb.jhu.edu/software/centrifuge/manual.shtml#database-download-and-index-building). 

- Centrifuge itself - Included in minotourapp/extras

- Set the environment variable MT_CENTRIFUGE_INDEX (the path to indexes without the suffix) to the location of your chosen Index and the centrifuge executable in the envs.sh file.

    export MT_CENTRIFUGE_INDEX="/path/to/centrifuge_indexes/p_compressed"

    The MT_CENTRIFUGE_INDEX above points to the following indexes:
    
    * p_compressed.1.cf
    
    * p_compressed.2.cf
    
    * p_compressed.3.cf
    
    * p_compressed.4.cf

- minoTour uses the ete3 package for Id'ing species from their taxonomic IDs. The installation of ete3's database during the first Metagenomics analysis will happen automatically, but requires internet access.

 - Alternatively, to force the download (this step just need to be executed once), type the following instructions into the python shell in the minoTour virtual environment:
```bash
    import ete3

    from ete3 import NCBITaxa

    n = NCBITaxa()
```

## Raise the Niceness index allowed

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
The centrifuge process is then launched with a niceness index of -10.
```bash
   <username> - nice -10
```

## Setting up Twitter
minoTour can tweet users when certain criteria are met - such as average coverage reaching nX, or minKNOW issuing warnings. To do this a twitter API account is needed. The steps for setting one up can be found [on twitters dev website]( https://developer.twitter.com/en/application/use-case),
and once the details have been issued to an account, the following environment variables need to be filled in [envs.sh](envs.sh). A twitter account is required.
MT_TWITCONSUMER_KEY='<Twitter consumer key>'
MT_TWITCONSUMER_SECRET='<Twitter consumer secret>'
MT_TWITTOKEN='<Twitter token>'
MT_TWITTOKEN_SECRET='<Twitter token secret>' 
    