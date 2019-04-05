from __future__ import absolute_import, unicode_literals

import os
import subprocess
import sys

import numpy as np
import pandas as pd
from celery import task
from celery.utils.log import get_task_logger
from django.conf import settings

from alignment.models import PafRoughCov, PafStore, PafSummaryCov
from jobs.models import JobMaster
from reads.models import FastqRead
from reference.models import ReferenceInfo, ReferenceLine

from web.utils import parse_md_cg_pafline, parseMDPAF, parse_mdpaf_alex

# Set up the logger to write to logging file
logger = get_task_logger(__name__)


def call_fetch_reads_alignment(runs, chunk_size, last_read):
    """
    Call the fetch reads function to create a fastq for the process
    :param runs: List of all runs on the flowcell
    :param chunk_size: The target number of reads we want
    :param last_read: The previous last read we took
    :return:
    """
    # a list of the fastqs object to pass into the mapping function
    fasta_objects = list()
    # Initialise the data frame
    fastq_df_barcode = pd.DataFrame()
    while True:
        # Call fetch_reads_alignment to actually query the database
        reads, last_read, read_count, fastasmchunk = fetch_reads_alignment(runs, chunk_size, last_read)
        # Add fasta_objects chunk to the list
        fasta_objects += fastasmchunk
        # Append the reads_df to the fastq_df
        fastq_df_barcode = fastq_df_barcode.append(reads)
        # If we have enough reads (more than chunk_size) or no reads
        if len(fastq_df_barcode) >= chunk_size or len(reads) == 0:
            break
    # Update the read count with the number of reads we just fetched
    read_count = len(fastq_df_barcode)
    return fastq_df_barcode, last_read, read_count, fasta_objects


def fetch_reads_alignment(runs, chunk_size, last_read):
    """
    Query the database for runs from a flowcell
    :param runs: The list of runs objects for this flowcell
    :param chunk_size: The target number of reads to have pulled back
    :param last_read: The id of the last read
    :return:
    """
    # Store the first read id and which run it's from this dict
    countsdict = dict()
    # The list of fastqread objects
    fasta_objects = list()
    # loop through the runs
    for run in runs:
        # fastqs = FastqRead.objects.filter(run=run, id__gt=int(last_read)).first()
        # Query the database for the reads objects
        fastqs = FastqRead.objects.values_list('id').filter(run=run, id__gt=int(last_read)).first()
        # If there are fastqs store the object
        if fastqs:
            # Store the first id you get and the urn it's from
            countsdict[fastqs[0]] = run
    # initialise the count
    count = 1
    # initialise the a dataframe to append reads to
    fastq_df_barcode = pd.DataFrame()
    # If you have reads from more than one run
    if len(countsdict) > 1:
        # Start with first read entry from this run
        for entry in sorted(countsdict):
            # If we haven't done more iterations than we have dict entries
            if count < len(countsdict):
                # Do the query
                fastqs = FastqRead.objects.filter(run=countsdict[entry],
                                                  id__gt=int(last_read),
                                                  id__lt=int(sorted(countsdict)[count]),
                                                  )[:chunk_size]
                # Append new results to barcode
                fastq_df_barcode = fastq_df_barcode.append(pd.DataFrame.from_records(
                    fastqs.values("read_id", "barcode_name", "sequence", "id")))

                # add to list of the fastQRead objects new objects
                fasta_objects += list(fastqs)
                # Set the new last read id
                last_read = fastq_df_barcode.id.max()
                # If we have enough reads, stop
                if len(fastq_df_barcode) >= chunk_size:
                    break
            count += 1
    # If we only have one run
    elif len(countsdict) == 1:
        """TODO This is risky and it breaks the logic - we may end up skipping reads from other runs"""
        # The key is the first fastqreads object retrieved
        mykey = list(countsdict.keys())[0]
        # Get the fastq objects from the database
        fastqs = FastqRead.objects.filter(run=countsdict[mykey],
                                          id__gt=int(last_read), )[:chunk_size]
        # append the results to the fastq dataframe
        fastq_df_barcode = fastq_df_barcode.append(pd.DataFrame.from_records(
            fastqs.values("read_id", "barcode_name", "sequence", "id")))
        # Create a list of the fastQRead objects we will use in this objects
        fasta_objects += list(fastqs)
        # get the new last read id
        last_read = fastq_df_barcode.id.max()
    # get the read count of the new reads pulled out to later update the job_master total
    read_count = len(fastq_df_barcode)
    # Return everything to be used above
    return fastq_df_barcode, last_read, read_count, fasta_objects


@task()
def run_minimap2_alignment_by_job_master(job_master_id):
    # TODO this is currently completely unused - To Delete? Plus run_task_minimap2.py
    """
    The code to run the minimap alignment from the django command line
    :param job_master_id: The id of the JobMaster for this alignment
    :return:
    """
    # Get the JobMaster Django Model object
    job_master = JobMaster.objects.get(pk=job_master_id)
    # If there is a reference
    if job_master.reference:
        run_minimap2_alignment(
            job_master.id,
        )
    else:
        logger.info('Flowcell id: {} - Job master {} has no reference.'.format(job_master.flowcell.id, job_master.id))


@task()
def run_minimap2_alignment(job_master_id):
    """
    :param job_master_id: The Primary key of the JobMaster object
    :return:
    """
    # Get the JobMaster
    job_master = JobMaster.objects.get(pk=job_master_id)
    # Get the flowcell
    flowcell = job_master.flowcell
    # Get all the runs attached to that Flowcell
    runs = flowcell.runs.all()
    # The chunk size of reads from this flowcell to fetch from the database
    chunk_size = 8000

    fastq_df_barcode, last_read, read_count, fasta_objects = call_fetch_reads_alignment(runs, chunk_size,
                                                                                        job_master.last_read)
    logger.info(fastq_df_barcode)
    logger.info('Flowcell id: {} - Running minimap2 on flowcell {}'.format(flowcell.id, flowcell.name))
    logger.info('Flowcell id: {} - job_master_id {}'.format(flowcell.id, job_master.id))
    logger.info('Flowcell id: {} - reference {}'.format(flowcell.id, job_master.reference.name))
    logger.info('Flowcell id: {} - last_read {}'.format(flowcell.id, job_master.last_read))
    logger.info('Flowcell id: {} - read_count {}'.format(flowcell.id, job_master.read_count))
    logger.info('Flowcell id: {} - number of reads found {}'.format(flowcell.id, read_count))
    # If we have pulled back reads, call fasta
    if read_count > 0:
        last_read = align_reads(fasta_objects, fastq_df_barcode, job_master.id)

    # Update the JobMaster with the metadata after finishing this iteration
    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = False
    job_master.last_read = last_read
    job_master.read_count = job_master.read_count + read_count
    job_master.save()


def align_reads(fastas, fasta_df, job_master_id):
    """
    Align reads to a reference. Requires a job_master.
    
    :param fastas: A list of FastqRead objects that contain the sequence we will be running alignment on
    :param fasta_df: A df of the Sequence, ReadID, ID and barcode_name
    :param job_master_id: The Primary key of the JobMaster model object
    :return:
    """

    reference_location = getattr(settings, "REFERENCE_LOCATION", None)
    minimap2 = getattr(settings, "MINIMAP2", None)

    job_master = JobMaster.objects.get(pk=job_master_id)

    last_read = job_master.last_read
    reference_info_id = job_master.reference.id
    flowcell = job_master.flowcell

    if not minimap2:
        logger.error('Can not find minimap2 executable - stopping task.')
        return

    if last_read is None:
        last_read = 0

    reference_info = ReferenceInfo.objects.get(pk=reference_info_id)

    chromdict = dict()

    chromosomes = reference_info.referencelines.all()

    for chromosome in chromosomes:
        chromdict[chromosome.line_name] = chromosome

    read = ''
    fastq_dict = dict()

    for fastq in fastas:
        # if not fastq.sequence:
        #    continue

        read = read + '>{} \r\n{}\r\n'.format(fastq.read_id, fastq.sequence)

        fastq_dict[fastq.read_id] = fastq

        last_read = fastq.id

    cmd = '{} -x map-ont -t 1 --secondary=no {} -'.format(
        minimap2,
        os.path.join(reference_location, reference_info.filename)
    )

    logger.info('Flowcell id: {} - Calling minimap - {}'.format(flowcell.id, cmd))

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            stdin=subprocess.PIPE, shell=True)

    (out, err) = proc.communicate(input=read.encode("utf-8"))

    status = proc.wait()

    logger.info('Flowcell id: {} - Finished minimap - {}'.format(flowcell.id, cmd))

    paf = out.decode("utf-8")

    pafdata = paf.splitlines()

    logger.info('Flowcell id: {} - Found {} paf records'.format(flowcell.id, len(pafdata)))

    if len(pafdata) > 0:

        paf_summary_cov_dict = dict()

        bulk_paf = []

        bulk_paf_rough = []

        for line in pafdata:

            line = line.strip('\n')

            record = line.split('\t')

            fastq_read = fastq_dict[record[0]]

            typeid = fastq_read.type.id

            is_pass = fastq_read.is_pass

            barcode_name = fastq_read.barcode_name

            newpaf = PafStore(
                job_master=job_master,
                read=fastq_read,
            )

            newpafstart = PafRoughCov(
                job_master=job_master,
                flowcell=flowcell,
                read_type=fastq_read.type,
                barcode_name=fastq_read.barcode_name,
                is_pass=fastq_read.is_pass
            )

            newpafend = PafRoughCov(
                job_master=job_master,
                flowcell=flowcell,
                read_type=fastq_read.type,
                barcode_name=fastq_read.barcode_name,
                is_pass=fastq_read.is_pass
            )

            newpaf.reference = reference_info

            newpaf.qsn = fastq_read.read_id  # models.CharField(max_length=256)#1	string	Query sequence name
            newpaf.qsl = int(record[1])  # models.IntegerField()#2	int	Query sequence length
            newpaf.qs = int(record[2])  # models.IntegerField()#3	int	Query start (0-based)
            newpaf.qe = int(record[3])  # models.IntegerField()#4	int	Query end (0-based)
            newpaf.rs = record[4]  # models.CharField(max_length=1)#5	char	Relative strand: "+" or "-"
            newpaf.tsn = chromdict[record[5]]  # models.CharField(max_length=256)#6	string	Target sequence name
            newpaf.tsl = int(record[6])  # models.IntegerField()#7	int	Target sequence length
            newpaf.ts = int(record[7])  # models.IntegerField()#8	int	Target start on original strand (0-based)
            newpaf.te = int(record[8])  # models.IntegerField()#9	int	Target end on original strand (0-based)
            newpaf.nrm = int(record[9])  # models.IntegerField()#10	int	Number of residue matches
            newpaf.abl = int(record[10])  # models.IntegerField()#11	int	Alignment block length
            newpaf.mq = int(record[11])  # models.IntegerField()#12	int	Mapping quality (0-255; 255 for missing)

            newpafstart.reference = reference_info
            newpafstart.chromosome = chromdict[record[5]]
            newpafstart.p = int(record[7])
            newpafstart.i = 1

            newpafend.reference = reference_info
            newpafend.chromosome = chromdict[record[5]]
            newpafend.p = int(record[8])
            newpafend.i = -1

            bulk_paf_rough.append(newpafstart)
            bulk_paf_rough.append(newpafend)

            bulk_paf.append(newpaf)

            if reference_info not in paf_summary_cov_dict:
                paf_summary_cov_dict[reference_info] = dict()

            if chromdict[record[5]] not in paf_summary_cov_dict[reference_info]:
                paf_summary_cov_dict[reference_info][chromdict[record[5]]] = dict()

            if barcode_name not in paf_summary_cov_dict[reference_info][chromdict[record[5]]]:
                paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name] = dict()

            if typeid not in paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name]:
                paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid] = dict()

            if 'read' not in paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]:
                paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]['read'] = set()
                paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]['length'] = 0

            paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]['read'].add(record[0])
            paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]['length'] += int(
                record[3]) - int(
                record[2]) + 1

        PafStore.objects.bulk_create(bulk_paf, batch_size=1000)
        PafRoughCov.objects.bulk_create(bulk_paf_rough, batch_size=1000)

        paf_store_list = PafStore.objects.filter(
            job_master=job_master
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

        paf_store_gb.reset_index().apply(lambda row: save_paf_store_summary(job_master.id, row), axis=1)

    return last_read


def save_paf_store_summary(job_master_id, row):
    barcode_name = row['read__barcode_name'][0]
    reference_line_name = row['tsn__line_name'][0]
    total_length = row['length']['sum']
    read_list = row['qsn']['unique']

    job_master = JobMaster.objects.get(pk=job_master_id)

    reference_line = ReferenceLine.objects.filter(reference=job_master.reference).filter(line_name=reference_line_name)[
        0]

    paf_summary_cov, created = PafSummaryCov.objects.get_or_create(
        job_master_id=job_master_id,
        barcode_name=barcode_name,
        reference_line_name=reference_line_name,
        reference_line_length=reference_line.chromosome_length
    )

    paf_summary_cov.total_length = total_length
    paf_summary_cov.read_count = len(read_list)

    paf_summary_cov.save()


def calculate_exepected_benefit_2dot0(flowcell_id, job_master_id):
    REFERENCE_LOCATION = getattr(settings, "REFERENCE_LOCATION", None)
    MINIMAP2 = getattr(settings, "MINIMAP2", None)
    CHUNK_SIZE = 2000

    job_master = JobMaster.objects.get(pk=job_master_id)

    last_read = job_master.last_read
    reference_info_id = job_master.reference.id
    flowcell = job_master.flowcell

    if last_read is None:
        last_read = 0

    fastqs = FastqRead.objects.filter(flowcell_id=flowcell_id, id__gt=int(last_read))[:CHUNK_SIZE]

    if not MINIMAP2:
        logger.error('Can not find minimap2 executable - stopping task.')
        print('Can not find minimap2 executable - stopping task.')
        return

    reference_info = ReferenceInfo.objects.get(pk=reference_info_id)

    chromdict = dict()

    chromosomes = reference_info.referencelines.all()

    for chromosome in chromosomes:
        chromdict[chromosome.line_name] = chromosome

    read = ''
    fastq_dict = dict()

    for fastq in fastqs:
        read = read + '>{} \r\n{}\r\n'.format(fastq.read_id, fastq.sequence)

        fastq_dict[fastq.read_id] = fastq

    cmd = '{} -x map-ont -t 1 --secondary=no -c --MD {} -'.format(
        MINIMAP2,
        os.path.join(REFERENCE_LOCATION, reference_info.filename)
    )

    logger.info('Flowcell id: {} - Calling minimap - {}'.format(flowcell.id, cmd))
    print('Flowcell id: {} - Calling minimap - {}'.format(flowcell.id, cmd))

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            stdin=subprocess.PIPE, shell=True)

    (out, err) = proc.communicate(input=read.encode("utf-8"))

    status = proc.wait()

    logger.info('Flowcell id: {} - Finished minimap - {}'.format(flowcell.id, cmd))
    print('Flowcell id: {} - Finished minimap - {}'.format(flowcell.id, cmd))

    paf = out.decode("utf-8")

    pafdata = paf.splitlines()

    logger.info('Flowcell id: {} - Found {} paf records'.format(flowcell.id, len(pafdata)))
    print('Flowcell id: {} - Found {} paf records'.format(flowcell.id, len(pafdata)))

    if len(pafdata) > 0:

        paf_summary_cov_dict = dict()

        bulk_paf = []

        bulk_paf_rough = []

        for line in pafdata:

            line = line.strip('\n')

            record = line.split('\t')

            fastq_read = fastq_dict[record[0]]

            typeid = fastq_read.type.id

            is_pass = fastq_read.is_pass

            barcode_name = fastq_read.barcode_name

            newpaf = PafStore(
                job_master=job_master,
                read=fastq_read,
            )

            newpafstart = PafRoughCov(
                job_master=job_master,
                flowcell=flowcell,
                read_type=fastq_read.type,
                barcode_name=fastq_read.barcode_name,
                is_pass=fastq_read.is_pass
            )

            newpafend = PafRoughCov(
                job_master=job_master,
                flowcell=flowcell,
                read_type=fastq_read.type,
                barcode_name=fastq_read.barcode_name,
                is_pass=fastq_read.is_pass
            )

            newpaf.reference = reference_info

            newpaf.qsn = fastq_read.read_id  # models.CharField(max_length=256)#1	string	Query sequence name
            newpaf.qsl = int(record[1])  # models.IntegerField()#2	int	Query sequence length
            newpaf.qs = int(record[2])  # models.IntegerField()#3	int	Query start (0-based)
            newpaf.qe = int(record[3])  # models.IntegerField()#4	int	Query end (0-based)
            newpaf.rs = record[4]  # models.CharField(max_length=1)#5	char	Relative strand: "+" or "-"
            newpaf.tsn = chromdict[record[5]]  # models.CharField(max_length=256)#6	string	Target sequence name
            newpaf.tsl = int(record[6])  # models.IntegerField()#7	int	Target sequence length
            newpaf.ts = int(record[7])  # models.IntegerField()#8	int	Target start on original strand (0-based)
            newpaf.te = int(record[8])  # models.IntegerField()#9	int	Target end on original strand (0-based)
            newpaf.nrm = int(record[9])  # models.IntegerField()#10	int	Number of residue matches
            newpaf.abl = int(record[10])  # models.IntegerField()#11	int	Alignment block length
            newpaf.mq = int(record[11])  # models.IntegerField()#12	int	Mapping quality (0-255; 255 for missing)

            newpafstart.reference = reference_info
            newpafstart.chromosome = chromdict[record[5]]
            newpafstart.p = int(record[7])
            newpafstart.i = 1

            newpafend.reference = reference_info
            newpafend.chromosome = chromdict[record[5]]
            newpafend.p = int(record[8])
            newpafend.i = -1

            bulk_paf_rough.append(newpafstart)
            bulk_paf_rough.append(newpafend)

            bulk_paf.append(newpaf)

            if reference_info not in paf_summary_cov_dict:
                paf_summary_cov_dict[reference_info] = dict()

            if chromdict[record[5]] not in paf_summary_cov_dict[reference_info]:
                paf_summary_cov_dict[reference_info][chromdict[record[5]]] = dict()

            if barcode_name not in paf_summary_cov_dict[reference_info][chromdict[record[5]]]:
                paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name] = dict()

            if typeid not in paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name]:
                paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid] = dict()

            if 'read' not in paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]:
                paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]['read'] = set()
                paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]['length'] = 0

            paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]['read'].add(record[0])
            paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]['length'] += int(
                record[3]) - int(
                record[2]) + 1

        PafStore.objects.bulk_create(bulk_paf, batch_size=1000)
        PafRoughCov.objects.bulk_create(bulk_paf_rough, batch_size=1000)

        paf_store_list = PafStore.objects.filter(
            job_master=job_master
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

        paf_store_gb.reset_index().apply(lambda row: save_paf_store_summary(job_master.id, row), axis=1)

    #
    # Now reprocess the pafdata and read the cigar for expected benefit.
    #

    # Create dictionaries if not passed in.
    # These dictionaries somehow need to persist from one iteration to the next.
    if last_read == 0:

        # Dictionary to hold matches, mismatches and coverage per base.
        # Structure referencedict[chromosome][match/mismatch/coverage][position]=[int]
        reference_dict = dict()

        # Dictionary to hold benefit. Structure benefitdict[chromosome][position]=[float]
        benefitdict = dict()

        # Dictionary to hold rolling benefit, mean benefit and mask.
        # This is badly structured! The Masks are stored as Bools with position.
        # The rollingbenefit is stored as a float with position.
        # The Means are floats but a single value per chromosome.
        # Structure rollingdict[chromosome][Forward/Reverse/ForMean/RevMean/ForMask/RevMask][optional position]=float,float,bool
        rollingdict = dict()

    else:

        raise NotImplementedError

    readcounter = 0
    mean = 0
    readlen = 0

    # These are parameters but might be tuned in future.
    error = 0.1
    prior_diff = 0.0001

    if len(pafdata) > 0:

        for line in pafdata:

            paf_record = parse_md_cg_pafline(line)

            if paf_record.chromosome not in reference_dict.keys():
                reference_dict[paf_record.chromosome] = {

                    'match': {

                    },
                    'mismatch': {

                    },
                    'coverage': {

                    },
                }

            # Get the specific arrays of matches, mismatches, start mapping position, end mapping position,
            # mapping orientation, reference(chromosome), reference lenght, readlength
            # These values are for a single read.
            # readlength is the total length of the read, not the aligned length.
            mismatcharray, matcharray, mapstart, mapend, maporientation, reference, referencelength, readlength = parse_mdpaf_alex(
                line)

            print(mismatcharray)
            print(matcharray)
            print(mapstart)
            print(mapend)
            print(maporientation)
            print(reference)
            print(referencelength)
            print(readlength)

            sys.exit()

            # print ("mismatcharray")
            # print(mismatcharray)
            # print(mismatcharray.dtype)
            # print ("Matcharray")
            # print(matcharray)
            # print(matcharray.dtype)

            # Calculate the length of the aligned portion we are looking at.
            readlen += len(mismatcharray)

            # Increment counter for reads processed on this loop.
            readcounter += 1

            # Calculate the mean of reads seen so far.
            mean = readlen / readcounter

            # If the chromosome we have seen is not in the reference dictionary,
            # we need to build it and also calculate the initial mask and rolling dictionary.
            if reference not in referencedict:
                referencedict, rollingdict = add_chromosome(
                    referencedict,
                    rollingdict,
                    benefitdict,
                    reference,
                    referencelength,
                    mean,
                    error,
                    prior_diff
                )

            # # Here we add the mismatch and match scores to the existing scores in the larger reference dictionary. We only update those positions which have changed for this read.
            # referencedict[reference]["mismatch"][mapstart:mapstart+len(mismatcharray)] += mismatcharray
            # referencedict[reference]["match"][mapstart:mapstart+len(matcharray)] += matcharray
            #
            # """
            # ### At this point we should have a dictionary for the reference with the count of matches and mismatches
            # ### Now we need to add this to any previously seen reads.
            # ### Suggest something like:
            # ### for reference in referencedict:
            # ###     fetch old coverage data
            # ###     add new coverage data to old coverage data
            # ###     now calculate getbenefit over the entire reference:
            # # Here we calculate the new benefit for this reference.
            # ### Somewhere we need to specify error and prior_diff - these are variables we might want to change.
            # """
            # benefitdict[reference] = getbenefit(
            #     referencedict[reference]["match"],
            #     referencedict[reference]["mismatch"],
            #     error=error, prior_diff=prior_diff)
            # A = benefitdict[reference]
            #
            # #mean is the mean read length at this point - to work for the rolling sum it has to be an integer!
            # rollingdict[reference]["Forward"] = rolling_sum([A], n=int(np.floor(mean)), pad=True)[0]
            # rollingdict[reference]["Reverse"] = rolling_sum([A[::-1]], n=int(np.floor(mean)), pad=True)[0][::-1]
            # # Recalculate the mean values for the reference we are looking at.
            # rollingdict[reference]["ForMean"] = rollingdict[reference]["Forward"].mean()
            # rollingdict[reference]["RevMean"] = rollingdict[reference]["Reverse"].mean()
            # # Update the coverage.
            # referencedict[reference]["coverage"] = np.sum(
            #     [referencedict[reference]["match"], referencedict[reference]["mismatch"]], axis=0)
            # # Update the mask for sequencing
            # rollingdict[reference]["ForMask"] = check_mask(rollingdict[reference]["Forward"],
            #                                                rollingdict[reference]["ForMean"])
            # rollingdict[reference]["RevMask"] = check_mask(rollingdict[reference]["Reverse"],
            #                                                rollingdict[reference]["RevMean"])

    print(reference_dict)

    last_read = fastq.id


def calculate_expected_benefit(flowcell_id, job_master_id):
    REFERENCE_LOCATION = getattr(settings, "REFERENCE_LOCATION", None)
    MINIMAP2 = getattr(settings, "MINIMAP2", None)
    CHUNK_SIZE = 2000

    job_master = JobMaster.objects.get(pk=job_master_id)

    last_read = job_master.last_read
    reference_info_id = job_master.reference.id
    flowcell = job_master.flowcell

    fastqs = FastqRead.objects.filter(flowcell_id=flowcell_id, id__gt=int(last_read))[:CHUNK_SIZE]

    if not MINIMAP2:
        logger.error('Can not find minimap2 executable - stopping task.')
        return

    if last_read is None:
        last_read = 0

    reference_info = ReferenceInfo.objects.get(pk=reference_info_id)

    chromdict = dict()

    chromosomes = reference_info.referencelines.all()

    for chromosome in chromosomes:
        chromdict[chromosome.line_name] = chromosome

    read = ''
    fastq_dict = dict()

    for fastq in fastqs:
        read = read + '>{} \r\n{}\r\n'.format(fastq.read_id, fastq.sequence)

        fastq_dict[fastq.read_id] = fastq

        last_read = fastq.id

    cmd = '{} -x map-ont -t 1 --secondary=no {} -'.format(
        MINIMAP2,
        os.path.join(REFERENCE_LOCATION, reference_info.filename)
    )

    logger.info('Flowcell id: {} - Calling minimap - {}'.format(flowcell.id, cmd))

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            stdin=subprocess.PIPE, shell=True)

    (out, err) = proc.communicate(input=read.encode("utf-8"))

    status = proc.wait()

    logger.info('Flowcell id: {} - Finished minimap - {}'.format(flowcell.id, cmd))

    paf = out.decode("utf-8")

    pafdata = paf.splitlines()

    logger.info('Flowcell id: {} - Found {} paf records'.format(flowcell.id, len(pafdata)))

    if len(pafdata) > 0:

        paf_summary_cov_dict = dict()

        bulk_paf = []

        bulk_paf_rough = []

        for line in pafdata:

            line = line.strip('\n')

            record = line.split('\t')

            fastq_read = fastq_dict[record[0]]

            typeid = fastq_read.type.id

            is_pass = fastq_read.is_pass

            barcode_name = fastq_read.barcode_name

            newpaf = PafStore(
                job_master=job_master,
                read=fastq_read,
            )

            newpafstart = PafRoughCov(
                job_master=job_master,
                flowcell=flowcell,
                read_type=fastq_read.type,
                barcode_name=fastq_read.barcode_name,
                is_pass=fastq_read.is_pass
            )

            newpafend = PafRoughCov(
                job_master=job_master,
                flowcell=flowcell,
                read_type=fastq_read.type,
                barcode_name=fastq_read.barcode_name,
                is_pass=fastq_read.is_pass
            )

            newpaf.reference = reference_info

            newpaf.qsn = fastq_read.read_id  # models.CharField(max_length=256)#1	string	Query sequence name
            newpaf.qsl = int(record[1])  # models.IntegerField()#2	int	Query sequence length
            newpaf.qs = int(record[2])  # models.IntegerField()#3	int	Query start (0-based)
            newpaf.qe = int(record[3])  # models.IntegerField()#4	int	Query end (0-based)
            newpaf.rs = record[4]  # models.CharField(max_length=1)#5	char	Relative strand: "+" or "-"
            newpaf.tsn = chromdict[record[5]]  # models.CharField(max_length=256)#6	string	Target sequence name
            newpaf.tsl = int(record[6])  # models.IntegerField()#7	int	Target sequence length
            newpaf.ts = int(record[7])  # models.IntegerField()#8	int	Target start on original strand (0-based)
            newpaf.te = int(record[8])  # models.IntegerField()#9	int	Target end on original strand (0-based)
            newpaf.nrm = int(record[9])  # models.IntegerField()#10	int	Number of residue matches
            newpaf.abl = int(record[10])  # models.IntegerField()#11	int	Alignment block length
            newpaf.mq = int(record[11])  # models.IntegerField()#12	int	Mapping quality (0-255; 255 for missing)

            newpafstart.reference = reference_info
            newpafstart.chromosome = chromdict[record[5]]
            newpafstart.p = int(record[7])
            newpafstart.i = 1

            newpafend.reference = reference_info
            newpafend.chromosome = chromdict[record[5]]
            newpafend.p = int(record[8])
            newpafend.i = -1

            bulk_paf_rough.append(newpafstart)
            bulk_paf_rough.append(newpafend)

            bulk_paf.append(newpaf)

            if reference_info not in paf_summary_cov_dict:
                paf_summary_cov_dict[reference_info] = dict()

            if chromdict[record[5]] not in paf_summary_cov_dict[reference_info]:
                paf_summary_cov_dict[reference_info][chromdict[record[5]]] = dict()

            if barcode_name not in paf_summary_cov_dict[reference_info][chromdict[record[5]]]:
                paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name] = dict()

            if typeid not in paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name]:
                paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid] = dict()

            if 'read' not in paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]:
                paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]['read'] = set()
                paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]['length'] = 0

            paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]['read'].add(record[0])
            paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]['length'] += int(
                record[3]) - int(
                record[2]) + 1

        PafStore.objects.bulk_create(bulk_paf, batch_size=1000)
        PafRoughCov.objects.bulk_create(bulk_paf_rough, batch_size=1000)

        paf_store_list = PafStore.objects.filter(
            job_master=job_master
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

        paf_store_gb.reset_index().apply(lambda row: save_paf_store_summary(job_master.id, row), axis=1)

    #
    # here start the expected benefit part
    #
    cmd = '{} -x map-ont -t 1 --secondary=no -c {} -'.format(
        MINIMAP2,
        os.path.join(REFERENCE_LOCATION, reference_info.filename)
    )

    logger.info('Flowcell id: {} - Calling minimap (CIGAR string) - {}'.format(flowcell.id, cmd))

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            stdin=subprocess.PIPE, shell=True)

    (out, err) = proc.communicate(input=read.encode("utf-8"))

    status = proc.wait()

    logger.info('Flowcell id: {} - Finished minimap (CIGAR string) - {}'.format(flowcell.id, cmd))

    paf = out.decode("utf-8")

    pafdata = paf.splitlines()

    if len(pafdata) > 0:

        for line in pafdata:
            pafline = parse_md_cg_pafline(line)

            reference_line = ReferenceLine.objects.filter(reference_info=reference_info).filter(
                line_name=pafline.chromosome)

            mismatch_hold_array = np.zeros(reference_line.chromosome_length)
            match_hold_array = np.zeros(reference_line.chromosome_length)

            mismatch_array, match_array, mapstart, runstart = parseMDPAF(line)
            mismatch_hold_array[mapstart:mapstart + len(mismatch_array)] += mismatch_array
            match_hold_array[mapstart:mapstart + len(match_array)] += match_array

            coverage_array = np.sum([match_hold_array, mismatch_hold_array], axis=0)

            filename_mismatch_hold = '{}_{}_mismatch_hold.bin'.format(flowcell.id, job_master.id)
            filename_match_hold = '{}_{}_match_hold.bin'.format(flowcell.id, job_master.id)
            filename_coverage = '{}_{}_coverage.bin'.format(flowcell.id, job_master.id)

            with open(os.path.join(REFERENCE_LOCATION, filename_mismatch_hold), 'wb+') as fh:
                fh.write(mismatch_hold_array.data)

            with open(os.path.join(REFERENCE_LOCATION, filename_match_hold), 'wb+') as fh:
                fh.write(match_hold_array.data)

            with open(os.path.join(REFERENCE_LOCATION, filename_coverage), 'wb+') as fh:
                fh.write(coverage_array.data)

    logger.info('Flowcell id: {} - Found {} paf records (CIGAR string)'.format(flowcell.id, len(pafdata)))

    return last_read
