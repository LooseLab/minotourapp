from __future__ import absolute_import, unicode_literals
import subprocess
from celery import task
from celery.utils.log import get_task_logger
from django.conf import settings
from alignment.models import PafRoughCov, PafStore
from readuntil.functions_EB import *
from jobs.models import JobMaster
import pandas as pd
import numpy as np
import pickle
from web.tasks_alignment import call_fetch_reads_alignment, save_paf_store_summary
from io import StringIO
from pathlib import Path
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
        calculate_expected_benefit_3dot0_final(
            job_master.id,
        )
    else:
        logger.info('Flowcell id: {} - Job master {} has no reference.'.format(job_master.flowcell.id, job_master.id))


def populate_priors_dict(reference_count_dict, priors_dict, fp, genotypes):
    """
    Populate the priors dictionary using the reference fasta file
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
    for desc, name, seq, qual in readfq(fp):
        # Call multi array results to create a structured array of zeros the lneght of the reference,
        #  with 9 fields, A,C,G,T,D,I,IC,M,U
        reference_count_dict[name] = multi_array_results(len(seq))
        # create a dictionary under the priors dict, keyed to the reference name
        priors_dict[name] = {}
        # add a priors key to this nested dictionary, then create and fill the priors for this reference,
        #  there's one prior for each reference position
        priors_dict[name]["priors"] = initialise_priors_rory(
            seq.upper(), genotypes
        )
        # add a posteriors key to this nested dictionary, then create and fill the posteriors array
        priors_dict[name]["posteriors"] = initialise_posteriors_mod(
            priors_dict[name]["priors"]
        )

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
    lam = 200
    m = 198
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
        print("model not recognized or implemented yet")

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
    fasta_df_barcode, last_read, read_count, fasta_objects = call_fetch_reads_alignment(runs, chunk_size,
                                                                                        task.last_read)

    if not minimap2:
        logger.error('Can not find minimap2 executable - stopping task.')
        return

    reference_info = task.reference
    # Chromosome dict contains all the chromsome names keyed to all the info we have about that chromosome
    chromdict = dict()
    # Get the chromosome we have from this reference
    chromosomes = reference_info.referencelines.all()
    # For each chromosome
    for chromosome in chromosomes:
        chromdict[chromosome.line_name] = chromosome

    # Create the series that contains the read_id and seqeunce as a correctly formatted fasta string
    fasta_df_barcode["fasta"] = ">" + fasta_df_barcode["read_id"] + "\n" + fasta_df_barcode["sequence"]
    # Create one string formatted correctly for fasta with input data
    fasta_data = "\n".join(list(fasta_df_barcode["fasta"]))
    # dict where we will be storing the fastq objects keyed to read ID
    fastq_dict = {fasta.read_id: fasta for fasta in fasta_objects}
    # create the command we are calling minimap with
    cmd = '{} -x map-ont -t 1 --secondary=no -c --MD {} -'.format(
        minimap2,
        Path(reference_location) / reference_info.filename
    )

    logger.info('Flowcell id: {} - Calling minimap EB - {}'.format(flowcell.id, cmd))
    # Setup minimap call
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            stdin=subprocess.PIPE, universal_newlines=True)
    # SubProcess communicate actually calls the function
    (out, err) = proc.communicate(input=fasta_data)

    logger.info('Flowcell id: {} - Finished minimap EB - {}'.format(flowcell.id, cmd))
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

        previously_known_chromosomes = set(ExpectedBenefitChromosomes.objects
                                              .filter(task=task).values_list("chromosome__line_name", flat=True))

        reference_path = Path(reference_location) / reference_info.filename

        if reference_info.filename.endswith(".gz"):
            with gzip.open(reference_path, "rt") as fp:
                priors_dict, reference_count_dict = populate_priors_dict(reference_count_dict, priors_dict, fp,
                                                                         genotypes)
        else:
            with open(reference_path, "r") as fp:
                priors_dict, reference_count_dict = populate_priors_dict(reference_count_dict, priors_dict, fp,
                                                                         genotypes)

        for i, record in enumerate(parse_PAF(StringIO(paf), fields=["qsn",
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
                                                                    "mq", "tags"])):
            # Get the read information from the fasta dict

            fastq_read = fastq_dict[record.qsn]

            chromosome_name = chromdict[record.tsn]

            chromosomes_seen_now.add(chromosome_name.line_name)

            newpaf = PafStore(
                job_master=task,
                read=fastq_read,
            )

            newpafstart = PafRoughCov(
                job_master=task,
                flowcell=flowcell,
                read_type=fastq_read.type,
                barcode_name=fastq_read.barcode_name,
                is_pass=fastq_read.is_pass
            )

            newpafend = PafRoughCov(
                job_master=task,
                flowcell=flowcell,
                read_type=fastq_read.type,
                barcode_name=fastq_read.barcode_name,
                is_pass=fastq_read.is_pass
            )

            newpaf.reference = reference_info

            newpaf.qsn = record.qsn  # models.CharField(max_length=256)#1	string	Query sequence name
            newpaf.qsl = record.qsl  # models.IntegerField()#2	int	Query sequence length
            newpaf.qs = record.qs  # models.IntegerField()#3	int	Query start (0-based)
            newpaf.qe = record.qe  # models.IntegerField()#4	int	Query end (0-based)
            newpaf.rs = record.rs  # models.CharField(max_length=1)#5	char	Relative strand: "+" or "-"
            newpaf.tsn = chromosome_name  # models.CharField(max_length=256)#6	Foreign key	Target sequence name
            newpaf.tsl = record.tsl  # models.IntegerField()#7	int	Target sequence length
            newpaf.ts = record.ts  # models.IntegerField()#8	int	Target start on original strand (0-based)
            newpaf.te = record.te  # models.IntegerField()#9	int	Target end on original strand (0-based)
            newpaf.nrm = record.nrm  # models.IntegerField()#10	int	Number of residue matches
            newpaf.abl = record.abl  # models.IntegerField()#11	int	Alignment block length
            newpaf.mq = record.mq  # models.IntegerField()#12	int	Mapping quality (0-255; 255 for missing)

            newpafstart.reference = reference_info
            newpafstart.chromosome = chromosome_name
            newpafstart.p = record.ts
            newpafstart.i = 1

            newpafend.reference = reference_info
            newpafend.chromosome = chromosome_name
            newpafend.p = record.te
            newpafend.i = -1

            bulk_paf_rough.append(newpafstart)
            bulk_paf_rough.append(newpafend)

            bulk_paf.append(newpaf)

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
                record.tags.get("cg", None), fastq_read.sequence, mapping_length
            )
            # Adds the mismatch count to the temporary dict, parsing the MD string
            d["M"] = np.fromiter(parse_MD(record.tags.get("MD", None)), np.uint16)

            # If the chromosome is not already a key in chrom dict, add it, with a 1 by 9 multidimensinal array,
            #  filled with zeros, the length of the reference, it should be there already but best to check
            if record.tsn not in reference_count_dict:
                reference_count_dict[record.tsn] = multi_array_results(record.tsl)
            # Loop through the dictionary of results, through each base
            for base in d:
                # add the counts for that base to the array that exists in the reference count dict for that base,
                # at the correct location
                reference_count_dict[record.tsn][base][
                mapping_start:mapping_end
                ] += d[base].astype(reference_count_dict[record.tsn][base].dtype)

            # TODO this needs to be reset, here we set the U to True to indicate that a count at this position
            # ToDO has been changed this iteration
            reference_count_dict[record.tsn]["U"][mapping_start:mapping_end] = True

            # For each chromosome in the reference
        iteration_counter = 0
        # Get how many chromsomes there are
        chromosome_total = len(reference_count_dict.keys()) - 1

        base_result_dir_path = make_results_directory(flowcell.id, task.id)

        for chrom_key in reference_count_dict:

            # Write out the coverage
            coverage_path = Path(f"{base_result_dir_path}/coverage_{chrom_key}_{flowcell.id}_{task.id}.dat")
            counts_path = Path(f"{base_result_dir_path}/counts_{chrom_key}_{flowcell.id}_{task.id}.dat")
            mask_path_forward = Path(f"{base_result_dir_path}/mask_forward_{chrom_key}_{flowcell.id}_{task.id}.dat")
            mask_path_reverse = Path(f"{base_result_dir_path}/mask_reverse_{chrom_key}_{flowcell.id}_{task.id}.dat")
            scores_forward_path = Path(f"{base_result_dir_path}/scores_forward_{chrom_key}_{flowcell.id}_{task.id}.dat")
            scores_reverse_path = Path(f"{base_result_dir_path}/scores_reverse_{chrom_key}_{flowcell.id}_{task.id}.dat")
            cost_forward_path = Path(f"{base_result_dir_path}/cost_forward_{chrom_key}_{flowcell.id}_{task.id}.dat")
            cost_reverse_path = Path(f"{base_result_dir_path}/cost_reverse_{chrom_key}_{flowcell.id}_{task.id}.dat")
            benefit_path = Path(f"{base_result_dir_path}/benefits_{chrom_key}_{flowcell.id}_{task.id}.dat")
            fixed_ben_forward_path = Path(f"{base_result_dir_path}/fixed_benefits_forward{chrom_key}"
                                          f"_{flowcell.id}_{task.id}.dat")
            fixed_ben_reverse_path = Path(f"{base_result_dir_path}/fixed_benefits_reverse{chrom_key}"
                                          f"_{flowcell.id}_{task.id}.dat")

            ############################################################################
            # ###################### read in the Counts ############################ #
            # ##################### Combine extant values ############################ #
            ############################################################################
            max_values = 0
            if counts_path.exists():
                # Read the old coverage counts into memory

                with open(counts_path, "rb") as fh:

                    old_counts_array = np.fromfile(fh, dtype=[("A", np.uint16), ("C", np.uint16),
                                                              ("G", np.uint16), ("T", np.uint16), ("D", np.uint16,),
                                                              ("I", np.uint16), ("IC", np.uint16),
                                                              ("M", np.uint16), ("U", np.bool)])
                    # Add the sum of each position counts to the old array, by each base type
                    for name in reference_count_dict[chrom_key].dtype.names:
                        # Don't update U as True + True = 2, which would ruin everything
                        if name is not "U":
                            reference_count_dict[chrom_key][name] += old_counts_array[name]
                            a = reference_count_dict[chrom_key][name].max()
                            if a > max_values:
                                max_values = a

            ############################################################################
            # ######################## Get coverage values ########################### #
            ############################################################################

            # Add the sum of each position counts
            coverage = np.sum(np.array(reference_count_dict[chrom_key][["A", "C", "G", "T"]].tolist(),
                                       dtype=np.int16), axis=1, dtype=np.uint16)

            ############################################################################

            print(f"iteration is {iteration_counter} of {chromosome_total}")

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
                    print("finished this chromosome")
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
                ] = update_posteriors_notts(posteriors_chunk, counts_chunk, phi)
                # update start size to index of the beginning of the next chunk
                start_slice += chunk_size

            # we have update all the chromosome in the while loop, so now we update the benefits
            priors_dict[chrom_key]["benefits"] = position_benefit_rory(
                priors_dict[chrom_key]["posteriors"], phi
            )
            # print(f"Benefits {priors_dict[chrom_key]['benefits']}")
            # print_a = priors_dict[chrom_key]['benefits']
            # print(f"Benefits zero size is {print_a[print_a==0].size}")
            # print(f"Benefits NaN size is {print_a[print_a==np.NaN].size} \n")

            # Calculate the fixed read lengths benefit - could use some clarification on how this works
            read_benefits_f, read_benefits_r = calculate_read_benefits_fixed_read_length_rory(
                priors_dict[chrom_key]["benefits"], lam, m
            )
            # print(f"Benefits fixed forward {read_benefits_f}")
            # print_a = read_benefits_f
            # print(f"Benefits fixed zero size is {print_a[print_a==0].size}")
            # print(f"Benefits fixed NaN size is {print_a[print_a==np.NaN].size} \n")
            #
            # print(f"Benefits fixed rev {read_benefits_r}")
            # print_a = read_benefits_r
            # print(f"Benefits fixed zero size is {print_a[print_a==0].size}")
            # print(f"Benefits fixed NaN size is {print_a[print_a==np.NaN].size} \n")
            # Calculate the m part of the read benefits of read lengths
            m_read_benefits_f, m_read_benefits_r = calculate_read_benefits_fixed_read_length_rory(
                priors_dict[chrom_key]["benefits"], m, m
            )
            # Costs of reading the read
            costs_f, costs_r = calculate_costs_fixed_read_length_rory(
                reference_length, rho, m, lam
            )
            # Scores are updated for whole reference
            scores_f, scores_r = calculate_scores_rory(
                read_benefits_f, read_benefits_r, m_read_benefits_f, m_read_benefits_r, costs_f, costs_r
            )
            # Accept F and accept R are the mask
            accept_f, accept_r = find_strategy_uniform(scores_f,
                                                       scores_r,
                                                       read_benefits_f,
                                                       read_benefits_r,
                                                       m_read_benefits_f,
                                                       m_read_benefits_r,
                                                       costs_f, costs_r,
                                                       rho, m, alpha)

            # ################################################################ #
            # ################# Write out coverage array ################ #
            ####################################################################
            # Add the max value for each of the arrays

            with open(coverage_path, "wb") as fh:
                fh.write(coverage.data)
            #del coverage
            # ################################################################ #
            # ################### Write out mask arrays ################## #
            ####################################################################
            # Forward ##########################################################
            with open(mask_path_forward, "wb") as fhm:
                fhm.write(accept_f.data)
            #del accept_f
            # Reverse ##########################################################
            with open(mask_path_reverse, "wb") as fhm:
                fhm.write(accept_r.data)
            #del accept_r
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
            #del priors_dict[chrom_key]["benefits"]
            ####################################################################

            # ################################################################ #
            # ################### Write out costs arrays ##################### #
            ####################################################################

            with open(cost_forward_path, "wb") as fhcost:
                fhcost.write(costs_f.data)
            #del costs_f
            # Reverse ##########################################################

            with open(cost_reverse_path, "wb") as fhcost:
                fhcost.write(costs_r.data)
            #del costs_r
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

        with open(Path("/home/rory/data/test_dict.pickle"), "wb") as fh:
            pickle.dump(priors_dict, fh, pickle.HIGHEST_PROTOCOL)

        PafStore.objects.bulk_create(bulk_paf, batch_size=1000)

        PafRoughCov.objects.bulk_create(bulk_paf_rough, batch_size=1000)

        paf_store_list = PafStore.objects.filter(
            job_master=task
        ).values(
            'job_master__id',
            'read__barcode_name',
            'tsn__line_name',
            'qsn',
            'qs',
            'qe'
        )

        paf_store_df = pd.DataFrame.from_records(paf_store_list)

        paf_store_df['length'] = paf_store_df['qe'] - paf_store_df['qs']

        paf_store_gb = paf_store_df.groupby(
            ['job_master__id', 'read__barcode_name', 'tsn__line_name']).agg(
            {'qsn': ['unique'], 'length': ['sum']})

        paf_store_gb.reset_index().apply(lambda row: save_paf_store_summary(task, row), axis=1)

        # this time we have seen these chromosomes for the first time, if any

        chromosomes_seen_now = chromosomes_seen_now - previously_known_chromosomes

        for chromosome_name in chromosomes_seen_now:
            eb_chromosome = ExpectedBenefitChromosomes(task=task, chromosome=chromdict[chromosome_name])
            eb_chromosome.save()

    task.running = False
    task.last_read = last_read
    task.read_count += read_count
    task.iteration_count += 1
    task.tempfile_name = "/tmp/benefit_results_dict.pickle"
    task.save()
