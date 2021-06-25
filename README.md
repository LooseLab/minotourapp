# minoTour

minoTour is a **work in progress**, therefore features may change, and any bug reports/feature requests/pull requests will be gratefully recieved.

minoTour is a web-based real-time laboratory information management system (LIMS) for Oxford Nanopore Technology (ONT) sequencers built using the Django framework. minoTour can be hosted offline and locally on a laptop, or online and accessed remotely.  It aims to enable monitoring of Nanopore sequencers in real time, allowing for analysis of performance and basic downstream data analysis. minoTour can monitor the activity of a sequencer in real time independent of analysing basecalled files providing a breakdown and analysis of live sequencing metrics via integration with ONT’s minKNOW API and parsing of sequence files as they are generated. These features enable remote control of a sequencer, as well as tracking flow cell QC and all messages received by MinKNOW during a sequencing run. The client for monitoring and uploading data, minFQ, is available via PyPI. Simple analyses, such as base calling summaries, alignment and metagenomics can be completed in real-time or after sequencing. minoTour is open source (https://www.github.com/LooseLab/minotourapp). 

An example instance can be found at http://137.44.59.170. Credentials to log in to a deactivated account (one which cannot upload data) are: 

    Username: demo-account
    Password: Welcome-to-minotour

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
    * [Adding references](#adding-references)   
    * [Metagenomics analyses](#metagenomics-analyses)
    * [Adding validation sets](#add-validation-sets-for-metagenomics)
    * [Raise allowed Niceness index](#raise-the-niceness-index-allowed)
    * [Setting up a twitter API account](#setting-up-twitter)
* [Uploading data to minoTour](#uploading-data-to-minotour)

## <a name="preconf"></a>Preconfiguration

Minotour uses MySQL as a database backend, Celery (asynchronous workers responsible for running server tasks), Redis (a database in memory similar to memcached, required by Celery as a message broker), Flower (optional - a useful Celery task monitor tool).

The following tools must be installed:

* [Redis](https://redis.io/download#installation) - Minotour uses **Redis** as a memory cache system and message broker for Celery. Follow the linked installation instructions.

* [MySQL Community edition](https://dev.mysql.com/downloads/) - Minotour requires a running instance of mySQL server. To set up mySQL community edition follow the instructions found in the [official installation documentation](https://dev.mysql.com/doc/mysql-getting-started/en/). Another good tutorial can be found [here](https://www.digitalocean.com/community/tutorials/how-to-install-mysql-on-ubuntu-16-04). Once the server is up and running, you can either use the root user created during set up, setting it in the below environment variables file.   
 Otherwise you can create a user using the mysql shell, whilst logged in as the root user or as an admin user:
    ```sql
    CREATE USER 'minotour'@'localhost' IDENTIFIED BY '<password>';
    CREATE DATABASE minotour;
    GRANT ALL PRIVILEGES ON 'minotour'.* TO 'minotour'@'localhost' 
    ```

* [Python 3](https://www.python.org) - Minotour uses **Python >=3.6**, so please make sure it is available on your system.

To upload data, it is also necessary to install the upload client, [minFQ](https://github.com/LooseLab/minotourcli), wherever you are uploading data from.
## <a name="setup"></a> Setting up the minoTour environment


### Set up minoTour
To set up a development environment, first clone or download the code from [our github page](https://github.com/LooseLab/minotourapp.git "Looselab's github page").
```bash
git clone https://github.com/LooseLab/minotourapp.git
```

Checkout the development branch:
    
    cd minotourapp
    git checkout develop

To use the python package mysqlclient, it is necessary to have two dependencies installed, libmysqlclinet-dev and python3-dev. These can be installed with the following command:

    sudo apt-get install libmysqlclient-dev
    sudo apt-get install python3-dev

To use weasyprint for pdf making we need to install it's dependencies
```bash
sudo apt-get install build-essential python3-dev python3-pip python3-setuptools python3-wheel python3-cffi libcairo2 libpango-1.0-0 libpangocairo-1.0-0 libgdk-pixbuf2.0-0 libffi-dev shared-mime-info
```

Now create a virtual environment for the project dependencies and install them:

    cd /path/to/minotourapp/
    python3 -m venv minotourenv
    source minotourenv/bin/activate
    pip install -r requirements.txt
    
### Setting up the environment variable file

Setting up environment variables - many of minoTour's config parameters are stored in the environment. There is a file to be configured,
[envs.sh](envs.sh). The fields required and what they represent are **explained in the file**.
For the secret key, we recommend generating your own, which can be done as follows:
```bash
$ python manage.py shell -c 'from django.core.management import utils; print(utils.get_random_secret_key())'
```
And stored under MT_SECRET_KEY in envs.sh.

### Setting the media root
The media root is found in minotourapp/settings.py. This is the directory that references and minimap2 
indexes are saved to. For example -
```python
MEDIA_ROOT = "/home/user/data/"
```

### Skipping the saving of sequence data

At the cost of not being able to go back and run analyses outside of real time against your upload of data,
 you can enable the skipping of saving the sequence data that is uploaded. This has the benefit of speeding up data upload substantially,
 and drastically reducing disk space usage by minoTour. This setting can be enabled in the environment variables file by setting:  
 
 ```bash
    MT_SKIP_SAVING_SEQUENCE=1
```

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
Once installed, minoTour has a [bash script](cripts/run_minotour_tmux.sh) that runs all the tmux commands to open minoTour in Tmux. *__The file may need configuring to provide all the correct paths__*, but with paths set should run without any issue.
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

* Start Celery::
    ```bash
        cd /path/to/minotour/code/ 
        source minotourenv/bin/activate
        source envs.sh
        celery -A minotourapp worker -l info -f ~/data/logs/celery.log --concurrency $WORKER_COUNT -Ofair
    ```
  
* Start Celery worker for minimap2::
    ```bash
        cd /path/to/minotour/code/ 
        source minotourenv/bin/activate
        source envs.sh
        celery -A minotourapp worker -l info -f $MT_LOG_FOLDER/celery.log -n minimap2 --concurrency $MT_CELERY_MINIMAP2_WORKER_COUNT -Ofair -Q minimap --pool threads
    ```
  
* Start Celery beat for periodic tasks::
    ```bash
        cd /path/to/minotour/code/ 
        source minotourenv/bin/activate
        source envs.sh
        celery beat -A minotourapp -l info
    ```

* Start Flower::
    ```bash
        cd /path/to/minotour/code/
        source minotourenv/bin/activate
        source envs.sh
        flower -A minotourapp --port=5555
    ```

* Start Minotour::
    ```bash
        cd /path/to/minotour/code
        source minotourenv/bin/activate
        source envs.sh
        python manage.py runserver 8100
    ```

* Time to test - if all went well, you should be able to access the web interface on http://localhost:8100.


## Final configurations

### Adding References
In order to perform alignments, it is necessary for the user to provide references to align to. There are currently two ways of doing this in minoTour,
using the GUI reference addition feature, or on the command line.

**1. Using the GUI reference addition:**
    It is possible to add references in the settings side bar, on the tab manage references. By clicking the plus icon on the page in the top right of the references table, a drag and drop area is popped up.
    However **the maximum file size for this method is 256Mb**, for larger references the other method must be used. If an identical reference has already been uploaded, minoTour will reject the upload, and will indicate which reference to use.
    References can be uploaded privately (only visible to the account) or publicly, visible to all users.
    
**2. Using the command line:**
    There is a command line command for adding larger references. **It requires an account api-key**, found in the sidebar under the profile tab. If running minoTour in tmux, there is a window pre-prepared for this, called zshelly. You can switch to it by typing 
```bash
CTRL+B, then 3
```
The command itself to add references is:
```bash
python manage.py add_references -k <account api key> /path/to/reference/file/or/directory/of/files
```
more help for the command can be received by typing
```bash
python manage.py add_references -h
```
If you are not using tmux (why not?) in a terminal shell run:
```bash
cd /path/to/minotour/code
source minotourenv/bin/activate
source envs.sh
python manage.py add_references -k <account api key> /path/to/reference/file/or/directory/of/files
``` 
   
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

## Add validation sets for metagenomics
Todo add this information
Remember that the name of the gff file must have an exact matching reference already uploaded.

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

## Uploading data to minoTour
To upload data, in both base-called and minKnow metrics format, we use a separate upload script, [minFQ](https://github.com/LooseLab/minotourcli "minFQ github page"). Example commands and how to upload data are found on the minFQ github.

 
    
