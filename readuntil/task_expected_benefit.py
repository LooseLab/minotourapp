from __future__ import absolute_import, unicode_literals
import subprocess
from celery import task
from celery.utils.log import get_task_logger
from django.conf import settings
from alignment.models import PafRoughCov, PafStore
from web.functions_EB import parse_PAF, parse_CIGAR, parse_MD
from jobs.models import JobMaster
import pandas as pd
import numpy as np
import pickle
import os
from web.tasks_alignment import call_fetch_reads_alignment, save_paf_store_summary
from io import StringIO

# Set up the logger to write to logging file
logger = get_task_logger(__name__)

@task()
def calculate_expected_benefit_3dot0_final(task_id):
    """
    Calculate the matches, mismatches and expected benefit banter
    :param task_id:
    :return:
    """
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
        os.path.join(reference_location, reference_info.filename)
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

    if paf:

        bulk_paf = []

        bulk_paf_rough = []
        # See if we have previously pickled results
        try:
            with open("/tmp/benefit_results.pickle", "rb") as data_store_fh:
                benefit_results_dict = pickle.load(data_store_fh)
        except FileNotFoundError as e:
            print("We don't have a file yet! Creating one now...")
            benefit_results_dict = {}
        #########################################################
        # ########### Parse paf update paf cov etc ############ #
        #########################################################
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
            fastq_read = fastq_dict[record[0]]

            chromosome_name = chromdict[record.tsn]

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
            newpaf.tsn = chromosome_name  # models.CharField(max_length=256)#6	string	Target sequence name
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
            mapping_start = record.ts
            mapping_end = record.te
            mapping_length = mapping_end-mapping_start
            # d contains our 7 arrays in theory
            # d IS NOT MULTIDIMENSIONAL
            q, d = parse_CIGAR(record.tags.get("cg", None), fastq_read.sequence, mapping_length)
            # Adds the mismatch count to the temporary dict
            d['M'] = parse_MD(record.tags.get("MD", None))
            # Results array
            # results_array = np.array([])
            # for key in d:
            #     # Put the results into multidimensional array,
            #     if results_array:
            #         np.vstack((results_array, d[key]))

            # If the chromosome is not already a key in chrom dict, add it, with a 1 by 8 multidimensinal array,
            #  filled with zeros, the length of the reference
            if chromosome_name not in benefit_results_dict:
                # SO THIS COULD BE MULTIDIMENSIONAL
                # loop through results and sum them at the correct position
                benefit_results_dict[chromosome_name] = dict()
                for key in d:
                    benefit_results_dict[chromosome_name][key] = np.zeros(reference_info.length)
                    benefit_results_dict[chromosome_name][key][mapping_start:mapping_end] += d[key]
            # if it is in the benefit results dict, yay! We already have an iteration
            else:
                for i, key in enumerate(d):
                    benefit_results_dict[chromosome_name][key][mapping_start:mapping_end] += d[key]

        with open("/tmp/benefit_results_dict.pickle", "wb") as fh_to_write:
            pickle.dump(benefit_results_dict, fh_to_write)

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

        task.running = False
        task.last_read = last_read
        task.read_count += chunk_size
        task.iteration_count += 1
        task.tempfile_name = "/tmp/benefit_results_dict.pickle"
        task.save()
