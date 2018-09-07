from __future__ import absolute_import, unicode_literals

import os
import subprocess

from celery import task
from celery.utils.log import get_task_logger
from django.conf import settings
from django.db.models import F

from alignment.models import PafRoughCov, PafStore
from devices.models import Flowcell
from jobs.models import JobMaster
from reads.models import FastqRead
from reference.models import ReferenceInfo

logger = get_task_logger(__name__)


@task
def run_minimap2_alignment_by_job_master(job_master_id):

    job_master = JobMaster.objects.get(pk=job_master_id)

    run_minimap2_alignment(
        job_master.flowcell.id,
        job_master.id,
        job_master.reference.id,
        job_master.last_read
    )

@task
def run_minimap2_alignment(flowcell_id, job_master_id, reference_info_id, last_read):

    logger.info("---> Task run_minimap2_alignment")
    logger.info("---> Parameters")
    logger.info("---> flowcell_id: {}".format(flowcell_id))
    logger.info("---> job_master_id: {}".format(job_master_id))
    logger.info("---> reference: {}".format(reference_info_id))
    logger.info("---> last_read: {}".format(last_read))

    if last_read is None:
        last_read = 0

    try:

        job_master = JobMaster.objects.get(pk=job_master_id)
        job_master.running = True
        job_master.save()

        data_folder = getattr(settings, "REFERENCELOCATION", None)

        reference_info = ReferenceInfo.objects.get(pk=reference_info_id)

        reference_info_path = os.path.join(data_folder, reference_info.minimap2_index_file_location)

        chromosome_dict = dict()

        chromosome_list = reference_info.referencelines.all()

        for chromosome in chromosome_list:

            chromosome_dict[chromosome.line_name] = chromosome

        fastq_dict = dict()
        fastq_type_dict = dict()
        fastq_barcode = dict()

        read_fastq = ''

        fastq_list = FastqRead.objects.filter(run__flowcell_id=flowcell_id, id__gt=int(last_read))[:1000]

        for fastq in fastq_list:

            read_fastq = read_fastq + '>{} \r\n{}\r\n'.format(fastq.read_id, fastq.fastqreadextra.sequence)
            fastq_dict[fastq.read_id] = fastq
            fastq_type_dict[fastq.read_id] = fastq.type
            fastq_barcode[fastq.read_id] = fastq.barcode
            last_read = fastq.id

        minimap2_executable = getattr(settings, "MINIMAP2", None)

        if not minimap2_executable:

            print('Can not find minimap2 executable - stopping task.')

        cmd = '{} -x map-ont -t 4 --secondary=no {} -'.format(minimap2_executable, reference_info_path)

        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                stdin=subprocess.PIPE, shell=True)

        (out, err) = proc.communicate(input=read_fastq.encode("utf-8"))

        paf = out.decode("utf-8")

        pafdata = paf.splitlines()

        flowcell = Flowcell.objects.get(pk=flowcell_id)

        resultstore = dict()

        bulk_paf = []
        bulk_paf_rough = []

        for line in pafdata:

            line = line.strip('\n')
            record = line.split('\t')

            fastq_read_id = fastq_dict[record[0]]
            fastq_read_type = fastq_type_dict[record[0]]
            barcode = fastq_barcode[record[0]]

            newpaf = PafStore(
                job_master=job_master,
                flowcell=flowcell,
                fastqread=fastq_read_id,
                reference_line=record[5]
            )

            newpaf.qsn = record[0]  # models.CharField(max_length=256)#1	string	Query sequence name
            newpaf.qsl = int(record[1])  # models.IntegerField()#2	int	Query sequence length
            newpaf.qs = int(record[2])  # models.IntegerField()#3	int	Query start (0-based)
            newpaf.qe = int(record[3])  # models.IntegerField()#4	int	Query end (0-based)
            newpaf.rs = record[4]  # models.CharField(max_length=1)#5	char	Relative strand: "+" or "-"
            newpaf.tsn = record[5]  # models.CharField(max_length=256)#6	string	Target sequence name
            # newpaf.tsn = chromdict[record[5]]  # models.CharField(max_length=256)#6	string	Target sequence name
            newpaf.tsl = int(record[6])  # models.IntegerField()#7	int	Target sequence length
            newpaf.ts = int(record[7])  # models.IntegerField()#8	int	Target start on original strand (0-based)
            newpaf.te = int(record[8])  # models.IntegerField()#9	int	Target end on original strand (0-based)
            newpaf.nrm = int(record[9])  # models.IntegerField()#10	int	Number of residue matches
            newpaf.abl = int(record[10])  # models.IntegerField()#11	int	Alignment block length
            newpaf.mq = int(record[11])  # models.IntegerField()#12	int	Mapping quality (0-255; 255 for missing)

            bulk_paf.append(newpaf)

            newpaf_start = PafRoughCov(
                flowcell=flowcell,
                read_type=fastq_read_type,
                read_type_name=fastq_read_type.name,
                barcode=barcode,
                barcode_name=barcode.name,
                chromosome=chromosome_dict[record[5]],
                p=int(record[7]),
                i=1
            )

            newpaf_end = PafRoughCov(
                flowcell=flowcell,
                read_type=fastq_read_type,
                read_type_name=fastq_read_type.name,
                barcode=barcode,
                barcode_name=barcode.name,
                chromosome=chromosome_dict[record[5]],
                p=int(record[8]),
                i=-1
            )

            bulk_paf_rough.append(newpaf_start)
            bulk_paf_rough.append(newpaf_end)

            # if reference_info not in resultstore:
            #     resultstore[reference_info] = dict()
            #
            # if chromosome_dict[record[5]] not in resultstore[reference_info]:
            #     resultstore[reference_info][chromosome_dict[record[5]]] = dict()
            #
            # if typeid not in resultstore[reference_info][chromosome_dict[record[5]]]:
            #     resultstore[reference_info][chromosome_dict[record[5]]][typeid] = dict()
            #
            # if 'read' not in resultstore[reference_info][chromosome_dict[record[5]]][typeid]:
            #     resultstore[reference_info][chromosome_dict[record[5]]][typeid]['read'] = set()
            #     resultstore[reference_info][chromosome_dict[record[5]]][typeid]['length'] = 0
            #
            # resultstore[reference_info][chromosome_dict[record[5]]][typeid]['read'].add(record[0])
            # resultstore[reference_info][chromosome_dict[record[5]]][typeid]['length'] += int(record[3]) - int(
            #     record[2]) + 1

        PafStore.objects.bulk_create(bulk_paf)
        PafRoughCov.objects.bulk_create(bulk_paf_rough)

        """
        for ref in resultstore:
            for ch in resultstore[ref]:
                for bc in resultstore[ref][ch]:
                    for ty in resultstore[ref][ch][bc]:
                        summarycov, created2 = PafSummaryCov.objects.update_or_create(
                            flowcell=flowcell,
                            read_type=ty,
                            barcodegroup=bc,
                            reference=ref,
                            chromosome=ch,
                        )
                        summarycov.read_count += len(resultstore[ref][ch][bc][ty]['read'])
                        summarycov.cumu_length += resultstore[ref][ch][bc][ty]['length']
                        summarycov.save()
        """

    except Exception as exception:  # Todo: This method of detecting an error will potentially lead to duplication of analysis
        print('An error occurred when running this task.')
        print(exception)
        JobMaster.objects.filter(pk=job_master_id).update(running=False)

    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = False
    job_master.last_read = last_read
    job_master.read_count = F('read_count') + len(fastq_list)
    job_master.save()

