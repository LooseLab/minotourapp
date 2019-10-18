from __future__ import absolute_import, unicode_literals

import os
import subprocess

import pandas as pd
from celery import task
from celery.utils.log import get_task_logger
from django.conf import settings

from alignment.models import PafRoughCov, PafStore, PafSummaryCov
from jobs.models import JobMaster
from reads.models import FastqRead
from reference.models import ReferenceInfo, ReferenceLine

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
    read_count = fastq_df_barcode.shape[0]
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

    fasta_df_barcode, last_read, read_count, fasta_objects = call_fetch_reads_alignment(runs, chunk_size,
                                                                                        job_master.last_read)

    logger.info('Flowcell id: {} - Running minimap2 on flowcell {}'.format(flowcell.id, flowcell.name))
    logger.info('Flowcell id: {} - job_master_id {}'.format(flowcell.id, job_master.id))
    logger.info('Flowcell id: {} - reference {}'.format(flowcell.id, job_master.reference.name))
    logger.info('Flowcell id: {} - last_read {}'.format(flowcell.id, job_master.last_read))
    logger.info('Flowcell id: {} - read_count {}'.format(flowcell.id, job_master.read_count))
    logger.info('Flowcell id: {} - number of reads found {}'.format(flowcell.id, read_count))
    # If we have pulled back reads, call fasta
    if read_count > 0:
        align_reads(fasta_objects, job_master.id, fasta_df_barcode)

    # Update the JobMaster with the metadata after finishing this iteration
    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = False
    job_master.last_read = last_read
    job_master.read_count = job_master.read_count + read_count
    job_master.save()


def align_reads(fastas, job_master_id, fasta_df):
    """
    Align reads to a reference. Requires a job_master.

    :param fastas: A list of FastqRead objects that contain the sequence we will be running alignment on
    :param job_master_id: The Primary key of the JobMaster model object
    :param fasta_df: A dataframe of fastq sequences
    :return:
    """

    # The on disk directory where we store the references
    reference_location = getattr(settings, "REFERENCE_LOCATION", None)
    # The location of the mimimap2 executable
    minimap2 = getattr(settings, "MINIMAP2", None)
    # The JobMaster object
    job_master = JobMaster.objects.get(pk=job_master_id)
    # Unpack metadata about this job
    last_read = job_master.last_read
    reference_info_id = job_master.reference.id
    flowcell = job_master.flowcell
    #
    if not minimap2:
        logger.error('Can not find minimap2 executable - stopping task.')
        return
    # We're starting fresh boys
    if last_read is None:
        last_read = 0
    # The information about the reference for this alignment
    reference_info = ReferenceInfo.objects.get(pk=reference_info_id)
    # A dictionary with all the results for each chromosome
    chromdict = dict()
    # Get all the chromosomes for this reference
    chromosomes = reference_info.referencelines.all()
    # For each Chromosome in the chromosomes, set the name of the chromosome as key, and the object as the value
    for chromosome in chromosomes:
        chromdict[chromosome.line_name] = chromosome

    # Create the series that contains the read_id and seqeunce as a correctly formatted fasta string
    fasta_df["fasta"] = ">" + fasta_df["read_id"] + "\n" + fasta_df["sequence"]
    # Create one string formatted correctly for fasta with input data
    fasta_data = "\n".join(list(fasta_df["fasta"]))

    # dict where we will be storing the fastq objects keyed to read ID
    fastq_dict = {fasta.read_id: fasta for fasta in fastas}
    # Create the minimap command
    cmd = '{} -x map-ont -t 1 --secondary=no {} -'.format(
        minimap2,
        os.path.join(reference_location, reference_info.filename)
    )

    logger.info('Flowcell id: {} - Calling minimap - {}'.format(flowcell.id, cmd))
    # Use subprocess to run minimap
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            stdin=subprocess.PIPE, universal_newlines=True)
    # Run the command and get a tuple of output
    (out, err) = proc.communicate(input=fasta_data)
    if err:
        logger.info(err)

    logger.info('Flowcell id: {} - Finished minimap - {}'.format(flowcell.id, cmd))
    # Get the results of the Paf output
    paf = out
    # Split the lines of the paf file
    pafdata = paf.splitlines()

    logger.info('Flowcell id: {} - Found {} paf records'.format(flowcell.id, len(pafdata)))
    # If we have paf data
    if len(pafdata) > 0:

        bulk_paf = []

        bulk_paf_rough = []
        # For each line in the paf data
        for line in pafdata:
            # Strip the new line char from the file
            line = line.strip('\n')
            # Split the line on the tab character
            record = line.split('\t')
            # Get the fastqread object to use in the Foreign key field below
            fastq_read = fastq_dict[record[0]]

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

            chromosome_line_name = chromdict[record[5]]
            # Set the fields on the PafStore
            newpaf.qsn = fastq_read.read_id  # models.CharField(max_length=256)#1	string	Query sequence name
            newpaf.qsl = int(record[1])  # models.IntegerField()#2	int	Query sequence length
            newpaf.qs = int(record[2])  # models.IntegerField()#3	int	Query start (0-based)
            newpaf.qe = int(record[3])  # models.IntegerField()#4	int	Query end (0-based)
            newpaf.rs = record[4]  # models.CharField(max_length=1)#5	char	Relative strand: "+" or "-"
            newpaf.tsn = chromosome_line_name  # models.CharField(max_length=256)#6	string	Target sequence name
            newpaf.tsl = int(record[6])  # models.IntegerField()#7	int	Target sequence length
            newpaf.ts = int(record[7])  # models.IntegerField()#8	int	Target start on original strand (0-based)
            newpaf.te = int(record[8])  # models.IntegerField()#9	int	Target end on original strand (0-based)
            newpaf.nrm = int(record[9])  # models.IntegerField()#10	int	Number of residue matches
            newpaf.abl = int(record[10])  # models.IntegerField()#11	int	Alignment block length
            newpaf.mq = int(record[11])  # models.IntegerField()#12	int	Mapping quality (0-255; 255 for missing)

            # Set the fields on the paf rough cov object
            newpafstart.reference = reference_info
            newpafstart.chromosome = chromosome_line_name
            newpafstart.p = int(record[7])
            newpafstart.i = 1

            # Set the end fields
            newpafend.reference = reference_info
            newpafend.chromosome = chromosome_line_name
            newpafend.p = int(record[8])
            newpafend.i = -1

            bulk_paf_rough.append(newpafstart)
            bulk_paf_rough.append(newpafend)

            bulk_paf.append(newpaf)

        PafStore.objects.bulk_create(bulk_paf, batch_size=1000)
        PafRoughCov.objects.bulk_create(bulk_paf_rough, batch_size=1000)

        paf_store_list = PafStore.objects.filter(job_master=job_master).values(
            'read__barcode_name',
            'tsn__line_name',
            'qsn',
            'qs',
            'qe'
        )
        # Create a DataFrame of the PafStore results from above, plus previous
        paf_store_df = pd.DataFrame.from_records(paf_store_list)
        # Minus the query start from the query end
        paf_store_df['length'] = paf_store_df['qe'] - paf_store_df['qs']
        # Do some maths?
        paf_store_gb = paf_store_df.groupby(
            ['read__barcode_name', 'tsn__line_name']).agg({'qsn': ['unique'], 'length': ['sum']})
        # Apply a saving function down the dataframe to save objects
        paf_store_gb.reset_index().apply(lambda row: save_paf_store_summary(job_master, row), axis=1)

    return last_read


def save_paf_store_summary(job_master, row):
    """
    Save the Paf Store Summaries
    :param job_master: The job master object
    :param row: The row of the dataframe
    :return:
    """

    barcode_name = row['read__barcode_name'][0]
    reference_line_name = row['tsn__line_name'][0]
    total_length = row['length']['sum']
    read_list = row['qsn']['unique']
    # The reference line for this object
    reference_line = ReferenceLine.objects.filter(reference=job_master.reference).filter(line_name=reference_line_name)[
        0]
    paf_summary_cov, created = PafSummaryCov.objects.get_or_create(
        job_master_id=job_master.id,
        barcode_name=barcode_name,
        reference_line_name=reference_line_name,
        reference_line_length=reference_line.chromosome_length
    )

    paf_summary_cov.total_length = total_length
    paf_summary_cov.read_count = len(read_list)

    paf_summary_cov.save()
