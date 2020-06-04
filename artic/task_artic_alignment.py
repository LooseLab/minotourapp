from __future__ import absolute_import, unicode_literals
from copy import deepcopy

import gzip
import os
import pickle
import subprocess
from collections import defaultdict
from io import StringIO

from celery import task
from celery.utils.log import get_task_logger

from alignment.models import PafSummaryCov
from alignment.tasks_alignment import call_fetch_reads_alignment
from minotourapp.utils import get_env_variable
from reads.models import JobMaster, JobType, Barcode, FastqReadType
from readuntil.functions_EB import *
from readuntil.models import ExpectedBenefitChromosomes

logger = get_task_logger(__name__)


# TODO split artic task into a sexy chain, rather than blocking a worker for so long

@task()
def run_artic_command(base_results_directory, barcode_name, job_master_pk):
    """
    Run the artic pipeline in this first pass method of dirty bash script
    Parameters
    ----------
    base_results_directory: str
        Path to the fastq file for this barcodes reads
    barcode_name: str
        The name of the barcode in string
    job_master_pk: int¯
        Primary key of the job master for this barcode
    Returns
    -------

    """
    jm = JobMaster.objects.get(pk=job_master_pk)
    jm.running = True
    jm.save()
    # Path to barcode fastq
    fastq_path = f"{base_results_directory}/{barcode_name}/{barcode_name}.fastq"

    # TODO get the barcode from the posix path for the sample name
    logger.info(fastq_path)
    scheme_dir = get_env_variable("MT_ARTIC_SCEHEME_DIR")
    os.chdir(f"{base_results_directory}/{barcode_name}")
    cmd = ["bash", "-c",
           f"source $CONDA_PREFIX/etc/profile.d/conda.sh && conda activate artic-ncov2019 && artic minion --medaka --normalise 200 --threads 4 --scheme-directory {scheme_dir} --read-file {fastq_path} nCoV-2019/V1 {barcode_name}"]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    out, err = proc.communicate()
    if not out and err:
        logger.info("We are done here. Artic out!")
        logger.info(out)
        logger.info(err)
    else:
        logger.info("¯\_(ツ)_/¯")
        logger.warning(str(out))
        logger.warning(err)
    if out and err:
        raise out

    jm = JobMaster.objects.get(pk=job_master_pk)
    jm.running = False
    jm.complete = True
    jm.save()

    # TODO delete the files we don't need here


def save_artic_job_masters(flowcell, barcode_name, reference_info):
    """
    Save the JobMasters that we need to manage the queue of tasks.
    Parameters
    ----------
    flowcell: reads.models.Flowcell
        Flowcell these runs are on
    barcode_name: str
        Barcode name for these reads
    reference_info: reference.models.ReferenceInfo
        Covid Reference
    Returns
    -------

    """

    job_type = JobType.objects.get(name="Run Artic")
    # TODO potential bug here where the barcode name on the reads is not the barcode name we save
    barcode_object = Barcode.objects.get(run__in=flowcell.runs.all(), name=barcode_name)
    job_master = JobMaster(job_type=job_type, reference=reference_info,
                           barcode=barcode_object, flowcell=flowcell)

    job_master.save()
    return

def make_barcoded_directories(barcodes, results_dir):
    """

    Parameters
    ----------
    barcodes: list of str
        List of barcode names
    results_dir: pathlib.PosixPath
        Path to the base results directory

    Returns
    -------
    None

    """
    for barcode in barcodes:
        path = results_dir / barcode
        if not path.exists():
            Path.mkdir(path)
            logger.info(f"Created directory {path} for barcode {barcode}.")
        else:
            logger.info(f"Directory {path} already exists.")


def add_barcode_to_tofire_file(barcode_name, directory):
    """
    Add a barcode to our file of barcodes we have fired the file for.
    Parameters
    ----------
    barcode_name: str
        Barcode name that we have enough coverage for to fire pipeline.
    directory: pathlib.PosixPath
        Absolute path to the results directory that we are storing this file in.
    Returns
    -------
    exit_code: int
        If 0 successfully added, else returns 1 and the error.
    """
    path = directory / "barcodes_to_fire.txt"

    with open(path, "a") as fh:
        try:
            fh.write(f"{barcode_name}\n")
            return 0
        except Exception as e:
            logger.error(str(e))
            return 1, e


def fetch_barcode_to_fire_list(directory):
    """

    Parameters
    ----------
    directory: pathlib.PosixPath
        the absolute path to the top level of th results directory of this task
    Returns
    -------
    barcode: list of str
        List of barcode to fire.

    """

    path = directory / "barcodes_to_fire.txt"
    try:
        with open(path, "r") as fh:
            barcodes = fh.read().split("\n")
    except FileNotFoundError as e:
        logger.info("File not found. First iteration? Continuing...")
        return []

    return (barcodes)


def make_results_directory_artic(flowcell_id, task_id, allow_create=True):
    """
    Make a results directory
    Parameters
    ----------
    flowcell_id: int
        Primary key of the flowcell entry in the database
    task_id: int
        Primary key of the task record in the database.
    allow_create: bool
        Allow the creation of the directory if it doesn't already exist

    Returns
    -------
    results_dir: pathlib.PosixPath
        PosixPath pointing to the results directory
    """
    environmental_results_directory = get_env_variable("MT_ARTIC_RESULTS_DIR")
    # environmental_results_directory = getattr(settings, "REFERENCE_LOCATION", None)

    artic_dir = Path(f"{environmental_results_directory}/artic/")
    if not artic_dir.exists() and allow_create:
        Path.mkdir(artic_dir)

    results_dir = Path(
        f"{environmental_results_directory}/artic/Temp_results"
    )

    if not results_dir.exists() and allow_create:
        Path.mkdir(results_dir)

    results_dir = Path(
        f"{environmental_results_directory}/artic/Temp_results/{flowcell_id}_{task_id}_artic"
    )

    if not results_dir.exists() and allow_create:
        Path.mkdir(results_dir)

    return results_dir


def populate_reference_count(
        reference_count_dict, fp, master_reference_count_path,
):
    """
    Populate the reference count dict so we can write it out with counts
    Parameters
    ----------
    reference_count_dict: dict
        Empty dictionary
    fp: file
        Filepath to the reference file
    master_reference_count_path: pathlib.PosixPath
        Path to the reference fasta on the server

    Returns
    -------
    reference_count_dict: dict
        The time of your life

    """
    # # Use the reference to populate the priors dict
    # If priors_dict and reference_count_dict for this dataset ALREADY EXISTS load them from archive

    if master_reference_count_path.exists():
        with open(master_reference_count_path, "rb") as f2:
            reference_count_dict = pickle.load(f2)
    else:
        for desc, name, seq, qual in readfq(fp):
            # Call multi array results to create a structured array of zeros the length of the reference,
            #  with 9 fields, A,C,G,T,D,I,IC,M,U
            reference_count_dict[name] = multi_array_results(len(seq))
            # create a dictionary under the priors dict, keyed to the reference name

        # Save reference_count_dict for next use
        with open(master_reference_count_path, "wb") as fh2:
            pickle.dump(reference_count_dict, fh2, pickle.HIGHEST_PROTOCOL)

    return reference_count_dict


def replicate_counts_array_barcoded(barcode_names, counts_dict):
    """
    Replicate the empty counts numpy array for each barcode
    Parameters
    ----------
    barcode_names: list of str
        the name of the barcode
    counts_dict: dict
        Dict keyed to chromosome/contig name. Value of the array for the counts.
        Array length of the number of bases in reference.

    Returns
    -------
    barcoded_counts_dict: dict
        Dictionary to store results. Keyed barcode name to a value of the reference count dictionary.
    None
    """
    # TODO this may be where the coverage ends up the same if the counts dict is shared, not unique to each barcode
    barcoded_counts_dict = {barcode: deepcopy(counts_dict) for barcode in barcode_names}
    return barcoded_counts_dict


@task()
def run_artic_pipeline(task_id, streamed_reads=None):
    """
    Run the artic pipeline on a flowcells amount of reads
    Parameters
    ----------
    task_id: int
        The pk of the task record in the database
    streamed_reads: list of dict
        List of reads in dictionary form. Default None.

    Returns
    -------

    """
    reference_location = getattr(settings, "REFERENCE_LOCATION", None)
    # The location of the mimimap2 executable
    minimap2 = getattr(settings, "MINIMAP2", None)

    task = JobMaster.objects.get(pk=task_id)

    if not task.reference:
        raise ValueError("Missing Reference file. Please sort out.")

    flowcell = task.flowcell

    avg_read_length = flowcell.average_read_length

    if avg_read_length == 0:
        logger.error(f"Average read length is zero Defaulting to 450, but this is an error.")
        avg_read_length = 450

    if not streamed_reads:
        runs = flowcell.runs.all()
        # The chunk size of reads from this flowcell to fetch from the database
        # aim for 100 megabases
        desired_yield = 100 * 1000000

        chunk_size = round(desired_yield / avg_read_length)

        logger.info(f"Fetching reads in chunks of {chunk_size} for alignment.")

        # Get the fastq objects
        fasta_df_barcode, last_read, read_count, fasta_objects = call_fetch_reads_alignment(
            runs, chunk_size, task.last_read, pass_only=True
        )

        if fasta_df_barcode.shape[0] == 0:
            logger.info("No fastq found. Skipping this iteration.")
            task.running = False
            task.save()
            return

    else:
        last_read = task.last_read
        fasta_df_barcode = pd.DataFrame(streamed_reads)
        fasta_df_barcode = fasta_df_barcode[fasta_df_barcode["is_pass"]]
        fasta_objects = streamed_reads

    read_count = fasta_df_barcode.shape[0]

    logger.info(f"Fetched reads.")
    logger.info(fasta_df_barcode.shape)

    min_read_length = 400

    max_read_length = 700

    fasta_df_barcode["sequence_length"] = fasta_df_barcode["sequence"].str.len()

    fasta_df_barcode = fasta_df_barcode[fasta_df_barcode["sequence_length"].between(
        min_read_length, max_read_length, inclusive=True)]

    logger.info(f"Reads after filtering for read lengths between"
          f" {min_read_length} - {max_read_length}: {fasta_df_barcode.shape[0]}")
    # TODO come back to look at this
    # fasta_df_barcode["barcode_name"] = np.where(fasta_df_barcode["barcode_name"].eq("No barcode"), "Unknown",
    #                                             fasta_df_barcode["barcode_name"])

    # ##Check if read_count is > 0 otherwise exit
    ## #We might want to check the minimum number of reads or the time since the last read was entered to do a smaller chunk?

    if read_count > 0:
        if not minimap2:
            logger.error("Can not find minimap2 executable - stopping task.")
            # This won't reset the tast?
            return

        reference_info = task.reference
        # Chromosome dict contains all the chromsome names keyed to all the info we have about that chromosome
        chromdict = dict()
        # Get the chromosome we have from this reference
        chromosomes = reference_info.referencelines.all()
        # For each chromosome
        for chromosome in chromosomes:
            chromdict[chromosome.line_name] = chromosome

        # Create the series that contains the read_id and sequence as a correctly formatted fasta string
        fasta_df_barcode["fasta"] = (
                ">" + fasta_df_barcode["read_id"] + "\n" + fasta_df_barcode["sequence"]
        )
        # Create one string formatted correctly for fasta with input data
        fasta_data = "\n".join(list(fasta_df_barcode["fasta"]))
        if streamed_reads:
            fastq_dict = {fasta.get("read_id"): fasta for fasta in fasta_objects}
        else:
            # dict where we will be storing the fastq objects keyed to read ID
            fastq_dict = {fasta.read_id: fasta for fasta in fasta_objects}
        # create the command we are calling minimap with
        cmd = "{} -x map-ont -t 1 --secondary=no -c --MD {} -".format(
            minimap2, Path(reference_location) / reference_info.filename
        )

        logger.info(
            "Flowcell id: {} - Calling minimap Artic - {}".format(flowcell.id, cmd)
        )
        # Setup minimap call
        proc = subprocess.Popen(
            cmd.split(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdin=subprocess.PIPE,
            universal_newlines=True,
        )
        # SubProcess communicate actually calls the function
        (out, err) = proc.communicate(input=fasta_data)

        if err:
            logger.error(err)

        logger.info(
            "Flowcell id: {} - Finished minimap EB - {}".format(flowcell.id, cmd)
        )
        # decode output from byte string
        paf = out
        logger.info("")

        if paf:
            reference_count_dict = {}

            chromosomes_seen_now = set()

            # Dictionary to store reads that have mapped under the correct barcode
            barcode_sorted_fastq_cache = defaultdict(list)

            reference_path = Path(reference_location) / reference_info.filename

            base_result_dir_path = make_results_directory_artic(flowcell.id, task.id)
            # We can save time by writing the empty numpy data struct and reloading it
            master_reference_count_path = Path(
                f"{base_result_dir_path}/{reference_info.name}_counts.pickle"
            )

            if reference_info.filename.endswith(".gz"):
                with gzip.open(reference_path, "rt") as fp:
                    reference_count_dict = populate_reference_count(
                        reference_count_dict, fp, master_reference_count_path
                    )
            else:
                with open(reference_path, "r") as fp:
                    reference_count_dict = populate_reference_count(
                        reference_count_dict, fp, master_reference_count_path
                    )

            barcodes = np.unique(fasta_df_barcode["barcode_name"])

            logger.info(f"Barcodes in this iteration are {barcodes}")

            barcoded_counts_dict = replicate_counts_array_barcoded(
                barcodes, reference_count_dict
            )

            paf_summary_cov_dict = {}

            make_barcoded_directories(barcodes, base_result_dir_path)

            logger.info(f"Parsing paf file. Please wait.")

            for i, record in enumerate(
                    parse_PAF(
                        StringIO(paf),
                        fields=[
                            "qsn",
                            "qsl",
                            "qs",
                            "qe",
                            "rs",
                            "tsn",
                            "tsl",
                            "ts",
                            "te",
                            "nrm",
                            "abl",
                            "mq",
                            "tags",
                        ],
                    )
            ):
                fastq_read = fastq_dict[record.qsn]

                if streamed_reads:
                    barcode_of_read = fastq_read.get("barcode_name", "Failed")
                    read_id = fastq_read.get("read_id")
                    sequence = fastq_read.get("sequence")
                    read_type = FastqReadType.objects.get(pk=fastq_read.get("type", 1))
                else:
                    barcode_of_read = fastq_read.barcode_name
                    read_id = fastq_read.read_id
                    sequence = fastq_read.sequence
                    read_type = fastq_read.type

                barcode_sorted_fastq_cache[barcode_of_read].append(f">{read_id}\n{sequence}")

                chromosome = chromdict[record.tsn]

                #########################################################
                # ########## Now for the parse cigar string ########### #
                #########################################################
                # Results array
                mapping_start = record.ts
                # The base at which the mapping ends
                mapping_end = record.te
                # The length of the mapping
                mapping_length = mapping_end - mapping_start

                if barcode_of_read not in paf_summary_cov_dict:
                    paf_summary_cov_dict[barcode_of_read] = {
                        "job_master": task,
                        "barcode_name": barcode_of_read,
                        "chromosome": chromosome,
                        "reference_line_length": chromosome.chromosome_length,
                        "yield": 0,
                        "read_count": 0,
                        "read_type": read_type
                    }

                paf_summary_cov = paf_summary_cov_dict[barcode_of_read]

                paf_summary_cov["yield"] += mapping_length

                paf_summary_cov["read_count"] += 1

                chromosomes_seen_now.add(chromosome.line_name)

                # Parse the cigar string, on the path file record, returns a dictionary of bases, with numpy arrays
                d = parse_cigar(
                    record.tags.get("cg", None), sequence, mapping_length,
                )
                # Adds the mismatch count to the temporary dict, parsing the MD string
                d["M"] = np.fromiter(parse_MD(record.tags.get("MD", None)), np.uint16)

                # Loop through the dictionary of results, through each base
                for base in d:
                    # add the counts for that base to the array that exists in the reference count dict for that base,
                    # at the correct location
                    barcoded_counts_dict[barcode_of_read][record.tsn][base][
                    mapping_start:mapping_end
                    ] += d[base].astype(reference_count_dict[record.tsn][base].dtype)

                # TODO this needs to be reset, here we set the U to True to indicate that a count at this position
                # ToDO has been changed this iteration, might offer speed up if we do this?
                # reference_count_dict[record.tsn]["U"][
                # mapping_start:mapping_end
                # ] = True

            # iterate the paf_summary_cov objects
            for paf_summary_barcode, paf_summary_cov in paf_summary_cov_dict.items():
                logger.info(f"Creating PafSummaryCov objects for {paf_summary_barcode}")
                paf_summary_cov_orm, created = PafSummaryCov.objects.get_or_create(
                    job_master=task,
                    barcode_name=paf_summary_cov["barcode_name"],
                    chromosome=paf_summary_cov["chromosome"],
                    reference_line_length=paf_summary_cov["reference_line_length"],
                    read_type=paf_summary_cov["read_type"],
                )

                paf_summary_cov_orm.total_length += paf_summary_cov["yield"]
                paf_summary_cov_orm.read_count += paf_summary_cov["read_count"]
                paf_summary_cov_orm.coverage = round(paf_summary_cov_orm.total_length
                                                                / paf_summary_cov_orm.reference_line_length, 2)
                paf_summary_cov_orm.average_read_length = round(paf_summary_cov_orm.total_length / paf_summary_cov_orm.read_count)

                paf_summary_cov_orm.save()

            # TODO only ever see one chromosome, so we can remove for loop?
            barcodes_already_fired = fetch_barcode_to_fire_list(base_result_dir_path)
            for chrom_key in chromosomes_seen_now:
                # Write out the coverage
                # TODO barcoding stuff goes here
                for barcode in barcodes:
                    if barcode in barcodes_already_fired:
                        continue
                    barcode_sorted_fastq_path = Path(f"{base_result_dir_path}/{barcode}/{barcode}.fastq")

                    coverage_path = Path(
                        f"{base_result_dir_path}/{barcode}/coverage_{chrom_key}_{flowcell.id}_{task.id}_{barcode}.dat"
                    )
                    counts_path = Path(
                        f"{base_result_dir_path}/{barcode}/counts_{chrom_key}_{flowcell.id}_{task.id}_{barcode}.dat"
                    )

                    ############################################################################
                    # ###################### read in the Counts ############################ #
                    # ##################### Combine extant values ############################ #
                    ############################################################################
                    max_values = 0
                    if counts_path.exists():
                        # Read the old coverage counts into memory

                        with open(counts_path, "rb") as fh:

                            old_counts_array = np.fromfile(
                                fh,
                                dtype=[
                                    ("A", np.uint16),
                                    ("C", np.uint16),
                                    ("G", np.uint16),
                                    ("T", np.uint16),
                                    ("D", np.uint16),
                                    ("I", np.uint16),
                                    ("IC", np.uint16),
                                    ("M", np.uint16),
                                    ("U", np.bool),
                                ],
                            )
                        # Add the sum of each position counts to the old array, by each base type
                        for name in barcoded_counts_dict[barcode][chrom_key].dtype.names:
                            # Don't update U as True + True = 2, which would ruin everything
                            if name is not "U":
                                barcoded_counts_dict[barcode][chrom_key][
                                    name
                                ] += old_counts_array[name]
                                a = barcoded_counts_dict[barcode][chrom_key][name].max()
                                if a > max_values:
                                    max_values = a

                    coverage = np.sum(
                        np.array(
                            barcoded_counts_dict[barcode][chrom_key][
                                ["A", "C", "G", "T"]
                            ].tolist(),
                            dtype=np.int16,
                        ),
                        axis=1,
                        dtype=np.uint16,
                    )

                    # ################################################################ #
                    # ################# Write out coverage array ################ #
                    ####################################################################
                    # Add the max value for each of the arrays

                    with open(coverage_path, "wb") as fh:
                        fh.write(coverage.data)

                    # ################################################################ #
                    # ################### Write out counts arrays ######################
                    ####################################################################

                    with open(counts_path, "wb") as fhc:
                        fhc.write(barcoded_counts_dict[barcode][chrom_key].data)

                    with open(barcode_sorted_fastq_path, "a") as fh:
                        fh.write("\n".join(barcode_sorted_fastq_cache[barcode]))

                    logger.info(
                        f"Avg Coverage at iteration {task.iteration_count} for barcode {barcode} is {coverage.mean()}"
                    )

                    number_bases_in_array = coverage.size
                    number_of_bases_250x = coverage[coverage >= 250].size
                    percent_250x = number_of_bases_250x / number_bases_in_array * 100
                    # 200x coverage
                    if np.all(coverage >= 200):
                        logger.info(
                            f"we would fire the code for the artic pipeline here for barcode {barcode}"
                        )
                        # run_artic_command(barcode_sorted_fastq_path, barcode)
                        add_barcode_to_tofire_file(barcode, base_result_dir_path)
                        save_artic_job_masters(flowcell, barcode, reference_info)

                        # task.running = False
                        # task.complete = True
                        # task.save()
                    elif percent_250x >= 95:

                        # So we have to set up the auto fire of the pipeline here if the task fails
                        logger.info(
                            f"we would fire the code for the artic pipeline here for barcode {barcode} due to 95% coverage reached"
                        )
                        # run_artic_command(barcode_sorted_fastq_path, barcode)
                        add_barcode_to_tofire_file(barcode, base_result_dir_path)
                        save_artic_job_masters(flowcell, barcode, reference_info)

                        # task.running = False
                        # task.complete = True
                        # task.save()
                        # Update last read id on task to start from on next iteration
                task.last_read = last_read
                task.iteration_count += 1

                logger.info("Finishing this batch of reads.")
                task.read_count += read_count
                task.running = False
                task.save()
