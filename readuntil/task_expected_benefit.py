from __future__ import absolute_import, unicode_literals

import gzip
import pickle
import subprocess
from io import StringIO
from pathlib import Path

import numpy as np
from celery import task
from celery.utils.log import get_task_logger

from alignment.tasks_alignment import call_fetch_reads_alignment
from reads.models import JobMaster
from readuntil.functions_EB import *
from readuntil.models import ExpectedBenefitChromosomes

# Set up the logger to write to logging file
logger = get_task_logger(__name__)


def run_eb_by_job_master(job_master_id):
    """
    The code to run the eb from the django command line
    :param job_master_id: The id of the JobMaster for this alignment
    :return:
    """
    # Get the JobMaster Django Model object
    job_master = JobMaster.objects.get(pk=job_master_id)
    # If there is a reference
    if job_master.reference:
        calculate_expected_benefit_3dot0_final(job_master.id)
    else:
        logger.info(
            "Flowcell id: {} - Job master {} has no reference.".format(
                job_master.flowcell.id, job_master.id
            )
        )


def populate_priors_dict(
    master_priors_path,
    master_reference_path,
    reference_count_dict,
    priors_dict,
    fp,
    genotypes,
):
    """
    Populate the priors dictionary using the reference fasta file
    :param master_priors_path: Path to a previously calculated set of priors for this reference.
    :param master_reference_path: Path to a previously calculated reference count for this reference.
    :param reference_count_dict: The counts dictionary for each base
    :type reference_count_dict: dict
    :param priors_dict: The dictionary that contains the priors for each chromosome in the reference
    :type priors_dict: dict
    :param fp: The file handle for the reference fasta file
    :type fp: file
    :param genotypes: The nucleotide bases we want to include - namely A, C, G, T and D
    :type genotypes: list
    :return: tuple of filled priors_dict and reference_count_dict
    """
    # # Use the reference to populate the priors dict
    ## If priors_dict and reference_count_dict for this dataset ALREADY EXISTS load them from archive
    ## This is sole

    if master_priors_path.exists():
        with open(master_priors_path, "rb") as f:
            priors_dict = pickle.load(f)
        with open(master_reference_path, "rb") as f2:
            reference_count_dict = pickle.load(f2)
    else:
        for desc, name, seq, qual in readfq(fp):
            # Call multi array results to create a structured array of zeros the length of the reference,
            #  with 9 fields, A,C,G,T,D,I,IC,M,U
            reference_count_dict[name] = multi_array_results(len(seq))
            # create a dictionary under the priors dict, keyed to the reference name
            priors_dict[name] = {}
            # add a priors key to this nested dictionary, then create and fill the priors for this reference,
            #  there's one prior for each reference position
            priors_dict[name]["priors"] = initialise_priors(
                seq.upper(), genotypes
            )
            # add a posteriors key to this nested dictionary, then create and fill the posteriors array
            priors_dict[name]["posteriors"] = initialise_posteriors_mod(
                priors_dict[name]["priors"]
            )

        ##Save priors_dict and reference_count_dict for next use

        with open(master_priors_path, "wb") as fh:
            pickle.dump(priors_dict, fh, pickle.HIGHEST_PROTOCOL)
        with open(master_reference_path, "wb") as fh2:
            pickle.dump(reference_count_dict, fh2, pickle.HIGHEST_PROTOCOL)
    return priors_dict, reference_count_dict


@task()
def calculate_expected_benefit_3dot0_final(task_id):
    """
    Calculate the matches, mismatches and expected benefit banter
    :param task_id:
    :return:
    """

    model = "DNAhaplo"
    errP = 0.05
    errD = 0.03
    errWeird = 0.05
    theta = 0.009
    r = 0.11
    lam = 4500
    m = 95
    rho = 0.2
    alpha = 0.5

    if model == "DNAhaplo":

        genotypes = ["A", "C", "G", "T", "D"]
        G = len(genotypes)
        indexG = {"A": 0, "C": 1, "G": 2, "T": 3, "D": 4}
        # probability that the correct base is read
        noErrP = 1.0 - (errP + errD)
        # phi are the probabilities of reading a base conditional on the actual genotype:
        # see equation (2) although here it's extended to include deletions.
        phi = []
        for i in range(G - 1):
            phi.append([])
            for j in range(G - 1):
                if i == j:
                    phi[i].append(noErrP)
                else:
                    phi[i].append(errP / (G - 2))
            phi[i].append(errD)
        phi.append([])
        for i in range(G - 1):
            phi[-1].append(errWeird / (G - 1))
        phi[-1].append(1.0 - errWeird)
    else:
        logger.info("model not recognized or implemented yet")

        # Change phi to np array rather than list
    phi = np.array(phi)
    # The on disk directory where we store the references
    reference_location = getattr(settings, "REFERENCE_LOCATION", None)
    # The location of the mimimap2 executable
    minimap2 = getattr(settings, "MINIMAP2", None)

    task = JobMaster.objects.get(pk=task_id)

    flowcell = task.flowcell

    runs = flowcell.runs.all()
    # The chunk size of reads from this flowcell to fetch from the database
    chunk_size = 8000
    # Get the fastq objects
    fasta_df_barcode, last_read, read_count, fasta_objects = call_fetch_reads_alignment(
        runs, chunk_size, task.last_read
    )

    ###Check if read_count is > 0 otherwise exit
    ###We might want to check the minimum number of reads or the time since the last read was entered to do a smaller chunk?

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
            ">"
            + fasta_df_barcode["read_id"]
            + "\n"
            + fasta_df_barcode["sequence"]
        )
        # Create one string formatted correctly for fasta with input data
        fasta_data = "\n".join(list(fasta_df_barcode["fasta"]))
        # dict where we will be storing the fastq objects keyed to read ID
        fastq_dict = {fasta.read_id: fasta for fasta in fasta_objects}
        # create the command we are calling minimap with
        cmd = "{} -x map-ont -t 1 --secondary=no -c --MD {} -".format(
            minimap2, Path(reference_location) / reference_info.filename
        )

        logger.info(
            "Flowcell id: {} - Calling minimap EB - {}".format(
                flowcell.id, cmd
            )
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
            "Flowcell id: {} - Finished minimap EB - {}".format(
                flowcell.id, cmd
            )
        )
        # decode output from byte string
        paf = out
        logger.info("")
        if paf:

            bulk_paf = []

            bulk_paf_rough = []
            # See if we have previously pickled results

            #########################################################
            # ########### Parse paf update paf cov etc ############ #
            #########################################################
            # Count of bases mapped to each chromosome at base resolution
            reference_count_dict = {}
            # Priors and posteriors and benefits for each chromosome will be stored in this dict
            priors_dict = {}
            # Mask dictionary - Forward and reverse mask will be stored in this dictionary

            chromosomes_seen_now = set()

            previously_known_chromosomes = set(
                ExpectedBenefitChromosomes.objects.filter(
                    task=task
                ).values_list("chromosome__line_name", flat=True)
            )

            reference_path = Path(reference_location) / reference_info.filename

            base_result_dir_path = make_results_directory(flowcell.id, task.id)

            current_working_directory = Path.cwd()

            master_priors_path = Path(
                f"{current_working_directory}/readuntil/Temp_results/priors_dict_{reference_info.filename}"
            )
            master_reference_path = Path(
                f"{current_working_directory}/readuntil/Temp_results/reference_dict_{reference_info.filename}"
            )

            #### Spit this into two functions. One will write the file, the second will read it.
            #### The code here writes the file.

            ####

            if reference_info.filename.endswith(".gz"):
                with gzip.open(reference_path, "rt") as fp:
                    priors_dict, reference_count_dict = populate_priors_dict(
                        master_priors_path,
                        master_reference_path,
                        reference_count_dict,
                        priors_dict,
                        fp,
                        genotypes,
                    )
            else:
                with open(reference_path, "r") as fp:
                    priors_dict, reference_count_dict = populate_priors_dict(
                        master_priors_path,
                        master_reference_path,
                        reference_count_dict,
                        priors_dict,
                        fp,
                        genotypes,
                    )

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
                # Get the read information from the fasta dict

                fastq_read = fastq_dict[record.qsn]

                chromosome_name = chromdict[record.tsn]

                chromosomes_seen_now.add(chromosome_name.line_name)

                #########################################################
                # ########## Now for the parse cigar string ########### #
                #########################################################
                # Results array
                mapping_start = record.ts
                # The base at which the mapping ends
                mapping_end = record.te
                # The length of the mapping
                mapping_length = mapping_end - mapping_start
                # Parse the cigar string, on the path file record, returns a dictionary of bases, with numpy arrays
                d = parse_cigar(
                    record.tags.get("cg", None),
                    fastq_read.sequence,
                    mapping_length,
                )
                # Adds the mismatch count to the temporary dict, parsing the MD string
                d["M"] = np.fromiter(
                    parse_MD(record.tags.get("MD", None)), np.uint16
                )

                # If the chromosome is not already a key in chrom dict, add it, with a 1 by 9 multidimensinal array,
                #  filled with zeros, the length of the reference, it should be there already but best to check

                # Here check if reference_count_dict[record.tsn] exists - if NOT then load that from the file and populate  priors_dict for this chromosome from the file.

                if record.tsn not in reference_count_dict:
                    reference_count_dict[record.tsn] = multi_array_results(
                        record.tsl
                    )
                # Loop through the dictionary of results, through each base
                for base in d:
                    # add the counts for that base to the array that exists in the reference count dict for that base,
                    # at the correct location
                    reference_count_dict[record.tsn][base][
                        mapping_start:mapping_end
                    ] += d[base].astype(
                        reference_count_dict[record.tsn][base].dtype
                    )

                # TODO this needs to be reset, here we set the U to True to indicate that a count at this position
                # ToDO has been changed this iteration
                reference_count_dict[record.tsn]["U"][
                    mapping_start:mapping_end
                ] = True

                # For each chromosome in the reference
            iteration_counter = 0
            # Get how many chromsomes there are
            chromosome_total = len(reference_count_dict.keys()) - 1

            ## Change code to only process data for chromosomes we have seen
            # for chrom_key in reference_count_dict:
            for chrom_key in chromosomes_seen_now:

                # Write out the coverage
                coverage_path = Path(
                    f"{base_result_dir_path}/coverage_{chrom_key}_{flowcell.id}_{task.id}.dat"
                )
                counts_path = Path(
                    f"{base_result_dir_path}/counts_{chrom_key}_{flowcell.id}_{task.id}.dat"
                )
                mask_path_forward = Path(
                    f"{base_result_dir_path}/mask_forward_{chrom_key}_{flowcell.id}_{task.id}.dat"
                )
                mask_path_reverse = Path(
                    f"{base_result_dir_path}/mask_reverse_{chrom_key}_{flowcell.id}_{task.id}.dat"
                )
                scores_forward_path = Path(
                    f"{base_result_dir_path}/scores_forward_{chrom_key}_{flowcell.id}_{task.id}.dat"
                )
                scores_reverse_path = Path(
                    f"{base_result_dir_path}/scores_reverse_{chrom_key}_{flowcell.id}_{task.id}.dat"
                )
                cost_forward_path = Path(
                    f"{base_result_dir_path}/cost_forward_{chrom_key}_{flowcell.id}_{task.id}.dat"
                )
                cost_reverse_path = Path(
                    f"{base_result_dir_path}/cost_reverse_{chrom_key}_{flowcell.id}_{task.id}.dat"
                )
                benefit_path = Path(
                    f"{base_result_dir_path}/benefits_{chrom_key}_{flowcell.id}_{task.id}.dat"
                )
                fixed_ben_forward_path = Path(
                    f"{base_result_dir_path}/fixed_benefits_forward{chrom_key}"
                    f"_{flowcell.id}_{task.id}.dat"
                )
                fixed_ben_reverse_path = Path(
                    f"{base_result_dir_path}/fixed_benefits_reverse{chrom_key}"
                    f"_{flowcell.id}_{task.id}.dat"
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
                        for name in reference_count_dict[
                            chrom_key
                        ].dtype.names:
                            # Don't update U as True + True = 2, which would ruin everything
                            if name is not "U":
                                reference_count_dict[chrom_key][
                                    name
                                ] += old_counts_array[name]
                                a = reference_count_dict[chrom_key][name].max()
                                if a > max_values:
                                    max_values = a

                ############################################################################
                # ######################## Get coverage values ########################### #
                ############################################################################

                # Add the sum of each position counts
                coverage = np.sum(
                    np.array(
                        reference_count_dict[chrom_key][
                            ["A", "C", "G", "T"]
                        ].tolist(),
                        dtype=np.int16,
                    ),
                    axis=1,
                    dtype=np.uint16,
                )

                ############################################################################

                logger.info(
                    f"iteration is {iteration_counter} of {chromosome_total}"
                )

                iteration_counter += 1

                # How many bases to cover for this chromosome in each while loop iteration
                chunk_size = 500000

                reference_length = reference_count_dict[chrom_key].size
                # How many bases in reference are left
                remainder = reference_length
                # start slice
                start_slice = 0
                # Bool to decide whether or not there are chunks to update the posteriors
                chunks_left = True

                while chunks_left:
                    # Remainder is the number of positions in the reference we haven't updated
                    remainder -= chunk_size
                    # if we have a remainder of 0 or less, we've covered the entire reference
                    if remainder <= 0:
                        logger.info("finished this chromosome")
                        # Set the end slice for this final iteration as the last element index
                        end_slice = start_slice + -remainder
                        # No more chunks left, so break out of loop after this iteration
                        chunks_left = False
                    # We have positions left after this iteration
                    else:
                        # Set the end slice plus how many bases we have in this iteration
                        end_slice = start_slice + chunk_size

                    # Get the posteriors for this chunk of the reference in a numpy array
                    posteriors_chunk = priors_dict[chrom_key]["posteriors"][
                        start_slice:end_slice
                    ]
                    # Get the counts of the number of each of the 5 bases that are
                    # mapped to the positions in this reference chunk
                    counts_chunk = reference_count_dict[chrom_key][
                        ["A", "C", "G", "T", "D"]
                    ][start_slice:end_slice]
                    # Update the posteriors in this chunk
                    priors_dict[chrom_key]["posteriors"][
                        start_slice:end_slice
                    ] = update_posteriors(
                        posteriors_chunk, counts_chunk, phi
                    )
                    # update start size to index of the beginning of the next chunk
                    start_slice += chunk_size

                # we have update all the chromosome in the while loop, so now we update the benefits
                priors_dict[chrom_key]["benefits"] = position_benefit(
                    priors_dict[chrom_key]["posteriors"], phi
                )
                # logger.info(f"Benefits {priors_dict[chrom_key]['benefits']}")
                # logger.info_a = priors_dict[chrom_key]['benefits']
                # logger.info(f"Benefits zero size is {logger.info_a[logger.info_a==0].size}")
                # logger.info(f"Benefits NaN size is {logger.info_a[logger.info_a==np.NaN].size} \n")

                # Calculate the fixed read lengths benefit - could use some clarification on how this works
                read_benefits_f, read_benefits_r = calculate_read_benefits_fixed_read_length(
                    priors_dict[chrom_key]["benefits"], lam, m
                )
                # logger.info(f"Benefits fixed forward {read_benefits_f}")
                # logger.info_a = read_benefits_f
                # logger.info(f"Benefits fixed zero size is {logger.info_a[logger.info_a==0].size}")
                # logger.info(f"Benefits fixed NaN size is {logger.info_a[logger.info_a==np.NaN].size} \n")
                #
                # logger.info(f"Benefits fixed rev {read_benefits_r}")
                # logger.info_a = read_benefits_r
                # logger.info(f"Benefits fixed zero size is {logger.info_a[logger.info_a==0].size}")
                # logger.info(f"Benefits fixed NaN size is {logger.info_a[logger.info_a==np.NaN].size} \n")
                # Calculate the m part of the read benefits of read lengths
                m_read_benefits_f, m_read_benefits_r = calculate_read_benefits_fixed_read_length(
                    priors_dict[chrom_key]["benefits"], m, m
                )
                # Costs of reading the read
                costs_f, costs_r = calculate_costs_fixed_read_length(
                    reference_length, rho, m, lam
                )
                # Scores are updated for whole reference
                scores_f, scores_r = calculate_scores(
                    read_benefits_f,
                    read_benefits_r,
                    m_read_benefits_f,
                    m_read_benefits_r,
                    costs_f,
                    costs_r,
                )
                # Accept F and accept R are the mask
                accept_f, accept_r = find_strategy_uniform(
                    scores_f,
                    scores_r,
                    read_benefits_f,
                    read_benefits_r,
                    m_read_benefits_f,
                    m_read_benefits_r,
                    costs_f,
                    costs_r,
                    rho,
                    m,
                    alpha,
                )

                # ################################################################ #
                # ################# Write out coverage array ################ #
                ####################################################################
                # Add the max value for each of the arrays

                with open(coverage_path, "wb") as fh:
                    fh.write(coverage.data)
                # del coverage
                # ################################################################ #
                # ################### Write out mask arrays ################## #
                ####################################################################
                # Forward ##########################################################
                with open(mask_path_forward, "wb") as fhm:
                    fhm.write(accept_f.data)
                # del accept_f
                # Reverse ##########################################################
                with open(mask_path_reverse, "wb") as fhm:
                    fhm.write(accept_r.data)
                # del accept_r
                ####################################################################

                # ################################################################ #
                # ################### Write out counts arrays ######################
                ####################################################################

                with open(counts_path, "wb") as fhc:
                    fhc.write(reference_count_dict[chrom_key].data)

                ####################################################################
                # ################### write out benefits array ################### #
                ####################################################################

                with open(benefit_path, "wb") as fhb:
                    fhb.write(priors_dict[chrom_key]["benefits"].data)
                # del priors_dict[chrom_key]["benefits"]
                ####################################################################

                # ################################################################ #
                # ################### Write out costs arrays ##################### #
                ####################################################################

                with open(cost_forward_path, "wb") as fhcost:
                    fhcost.write(costs_f.data)
                # del costs_f
                # Reverse ##########################################################

                with open(cost_reverse_path, "wb") as fhcost:
                    fhcost.write(costs_r.data)
                # del costs_r
                # ################################################################ #
                # ################### Write out scores arrays #################### #
                ####################################################################

                with open(scores_forward_path, "wb") as fhscores:
                    fhscores.write(scores_f.data)

                # Reverse ##########################################################

                with open(scores_reverse_path, "wb") as fhscores:
                    fhscores.write(scores_r.data)

                # ################################################################ #
                # ################### Write out fixed benefit #################### #
                ####################################################################

                with open(fixed_ben_forward_path, "wb") as fhbf:
                    fhbf.write(read_benefits_f.data)

                # Reverse ##########################################################

                with open(fixed_ben_reverse_path, "wb") as fhbr:
                    fhbr.write(read_benefits_r.data)

            # with open(Path("/home/rory/data/test_dict.pickle"), "wb") as fh:
            #    pickle.dump(priors_dict, fh, pickle.HIGHEST_PROTOCOL)

            # this time we have seen these chromosomes for the first time, if any

            chromosomes_seen_now = (
                chromosomes_seen_now - previously_known_chromosomes
            )
            task = JobMaster.objects.get(pk=task_id)
            for chromosome_name in chromosomes_seen_now:
                eb_chromosome = ExpectedBenefitChromosomes(
                    task=task, chromosome=chromdict[chromosome_name]
                )
                eb_chromosome.save()

    task = JobMaster.objects.get(pk=task_id)
    task.running = False
    task.last_read = last_read
    task.read_count += read_count
    task.iteration_count += 1
    task.tempfile_name = "/tmp/benefit_results_dict.pickle"
    task.save()
