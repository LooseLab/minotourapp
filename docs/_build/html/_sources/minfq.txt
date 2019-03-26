###############
Minotour client
###############

**The client is dependent on python3.5 and above.** We recommend creating a virtual environment to contain the environment and avoid polluting your global environment.::

    python3 -m venv /path/to/envs/minfq

Activate the environment::

    source /path/to/envs/minfq/bin/activate

Upgrade pip to its latest version::

    pip install --upgrade pip

The client is available on PyPi::

    pip install minFQ

A development version is also available. This version is usually more up to date. Clone the client repository::

    git clone https://github.com/LooseLab/minotourcli.git

Or download and unzip the code found in the github `repository <https://github.com/LooseLab/minotourcli>`_.

To install the development version into the python virtual environment::

    pip install -e .

Check the client is installed into the virtual environment (make sure the environment is activated, using the above source command)::

    minFQ -h

You should see a helpful help page.

To upload data to the minotour application, you will need the API key of the registered user you wish the data to be kept under. This can be found in the profile section of this user.

To access the profile page, login using your username and password, click the username in the top nav bar, click the profile option. Copy the API key to clipboard.

An example minFQ command would be::

    minFQ -w /path/to/directory/containing_fastq/ -n <flowcell_name> -k <Api_key> -hn <Server address for minotour> -ip <minKNOW address> -p <Port number>

For docker the port is 10000. The host server address for minotour is localhost, or 127.0.0.1. The minKNOW ip address is usually 127.0.0.1. The flowcell name can be what you wish, and the watch directory is where the basecalled Fastq will be appearing.

minFQ scans recursively, so setting -w a top level directory for the flowcell data will find all the runs and Fastqs inside that directory.

For the development environment, this will be similar, but the port can vary.

For the minotour web server hosted by Nottingham, the host server address for minotour will be minotour.nottingham.ac.uk.

-----------
Config file
-----------

minFQ can be preconfigured with a lot of the options that stay the same using a config file, that must be present in the Current Working Directory that minFQ is being called from.

An example config file could be called minfq-posix.config, and any option can be configured using the -- name of that argument. --names are viewable by running::

    minFQ -h

.Config file syntax allows: key=value, flag=true, stuff=[a,b,c].

The contents of a config file may look like::

    key=b410a0c9729d92ac989509c695cc3ee66a749ec6
    port=10000
    hostname=127.0.0.1
    ip-address=127.0.0.1


