from __future__ import absolute_import, unicode_literals

import os
import subprocess

import pandas as pd
from celery import task
from celery.utils.log import get_task_logger
from django.conf import settings

from alignment.models import PafRoughCov, PafStore, PafSummaryCov
from jobs.models import JobMaster
from reads.models import FastqRead, Flowcell
from reference.models import ReferenceInfo, ReferenceLine

logger = get_task_logger(__name__)


@task
def run_minimap2_alignment_by_job_master(job_master_id):

    job_master = JobMaster.objects.get(pk=job_master_id)

    if job_master.reference:

        run_minimap2_alignment(
            job_master.flowcell.id,
            job_master.id,
            job_master.reference.id,
            job_master.last_read
        )

    logger.info('Flowcell id: {} - Job master {} has no reference.'.format(job_master.flowcell.id, job_master.id))


@task
def run_minimap2_alignment(flowcell_id, job_master_id, reference_info_id, last_read):

    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = True
    job_master.save()

    flowcell = Flowcell.objects.get(pk=flowcell_id)

    fastqs = FastqRead.objects.filter(run__flowcell_id=flowcell_id, id__gt=int(last_read))[:1000]

    logger.info('Flowcell id: {} - Running minimap2 on flowcell {}'.format(flowcell.id, flowcell.name))
    logger.info('Flowcell id: {} - job_master_id {}'.format(flowcell.id, job_master.id))
    logger.info('Flowcell id: {} - reference {}'.format(flowcell.id, job_master.reference.name))
    logger.info('Flowcell id: {} - last_read {}'.format(flowcell.id, job_master.last_read))
    logger.info('Flowcell id: {} - read_count {}'.format(flowcell.id, job_master.read_count))
    logger.info('Flowcell id: {} - number of reads found {}'.format(flowcell.id, len(fastqs)))

    if len(fastqs) > 0:

        REFERENCE_LOCATION = getattr(settings, "REFERENCE_LOCATION", None)
        MINIMAP2 = getattr(settings, "MINIMAP2", None)

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

        minimap2 = reference_info.minimap2_index_file_location

        read = ''
        fastq_dict = dict()
        fastqtypedict = dict()
        fastq_read_ispass_dict = dict()
        fastq_read_barcode_dict = dict()

        fastqbarcode=dict()

        for fastq in fastqs:

            read = read + '>{} \r\n{}\r\n'.format(fastq.read_id, fastq.fastqreadextra.sequence)

            fastq_dict[fastq.read_id] = fastq
            fastqtypedict[fastq.read_id] = fastq.type
            fastq_read_ispass_dict[fastq.read_id] = fastq.is_pass
            fastq_read_barcode_dict[fastq.read_id] = fastq.barcode_name
            fastqbarcode[fastq.read_id] = fastq.barcode.name

            last_read = fastq.id

        cmd = '{} -x map-ont -t 4 --secondary=no {} -'.format(
            MINIMAP2,
            os.path.join(REFERENCE_LOCATION, reference_info.filename)
        )

        logger.info('Flowcell id: {} - Calling minimap - {}'.format(flowcell.id, cmd))

        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                stdin=subprocess.PIPE, shell=True)

        (out, err) = proc.communicate(input=read.encode("utf-8"))

        status = proc.wait()

        paf = out.decode("utf-8")

        pafdata = paf.splitlines()

        paf_summary_cov_dict = dict()

        ### OK - we need to fix this for getting the group barcodes and not the individual barcodes.

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
            paf_summary_cov_dict[reference_info][chromdict[record[5]]][barcode_name][typeid]['length'] += int(record[3]) - int(
                record[2]) + 1

        PafStore.objects.bulk_create(bulk_paf)
        PafRoughCov.objects.bulk_create(bulk_paf_rough)

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

    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = False
    job_master.last_read = last_read
    job_master.read_count = job_master.read_count + len(fastqs)
    job_master.save()


def save_paf_store_summary(job_master_id, row):

    barcode_name = row['read__barcode_name'][0]
    reference_line_name = row['tsn__line_name'][0]
    total_length = row['length']['sum']
    read_list = row['qsn']['unique']

    job_master = JobMaster.objects.get(pk=job_master_id)

    reference_line = ReferenceLine.objects.filter(reference=job_master.reference).filter(line_name=reference_line_name)[0]

    paf_summary_cov, created = PafSummaryCov.objects.get_or_create(
        job_master_id=job_master_id,
        barcode_name=barcode_name,
        reference_line_name=reference_line_name,
        reference_line_length=reference_line.chromosome_length
    )

    paf_summary_cov.total_length = total_length
    paf_summary_cov.read_count = len(read_list)

    paf_summary_cov.save()
