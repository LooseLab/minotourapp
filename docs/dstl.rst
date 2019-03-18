#####################
Threat agent analysis
#####################

To start a target analysis from the client it is first necessary to add the references to the minotour web app. This is covered whilst setting up docker here - :doc:`docker`

The -j flags and -ts flags on the minotour client, minFQ, can be used to automatically start a job upon flowcell creation.

To see which jobs, references and validation sets are available to you, after setup of the client :doc: minfq ::

    source myenv/bin/activate
    minFQ --list.

We recommend using a config file, as covered in the client setup instructions :doc: minfq.

In this config file, set::

    -j Metagenomics
    -ts <Target_set_name>

Run the client upload, minFQ -n <Flowcell_name> -w /path/to/fastq/directory.

--------------------------------------
Start the analysis in the minotour app
--------------------------------------

It is possible to start the task in the minotour app, by navigating into the view of the desired flowcell, switching to the tabs task and selecting the Metagenomics task.

If a target set is not selected from the second drop down, no mapping will take place, and only Centrifuge based results will be shown.