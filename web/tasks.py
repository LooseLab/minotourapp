from __future__ import absolute_import, unicode_literals

import json
import os
import subprocess
import tempfile
from datetime import datetime, timedelta

from django.db.models import F

import pytz
import redis
from celery import task
from celery.utils.log import get_task_logger
from django.conf import settings
from django.contrib.auth.models import User
from django.core.cache import cache
from django.core.mail import send_mail
from django.core.paginator import Paginator
from django.db.models import Q
from django_mailgun import MailgunAPIError
from twitter import *

from alignment.models import PafStore,PafStore_transcriptome
from alignment.models import PafSummaryCov,PafSummaryCov_transcriptome
from alignment.models import SamStore
from communication.models import Message
from minikraken.models import MiniKraken, ParsedKraken
from reads.models import Barcode
from reads.models import ChannelSummary
from reads.models import FastqRead
from reads.models import HistogramSummary
from reads.models import JobMaster
from reads.models import JobType
from reads.models import MinIONRun, FastqReadType
from reads.models import RunStatisticBarcode
from reads.models import RunSummaryBarcode
from reference.models import ReferenceInfo

logger = get_task_logger(__name__)


def utcnow():
    return datetime.now(tz=pytz.utc)


def findbin(x, bin_width=900):
    return (x - x % bin_width) / bin_width


@task()
def run_monitor():

    logger.info('Running run_monitor celery task.')

    print("Rapid Monitor Called")

    minion_runs = MinIONRun.objects.filter(active=True).distinct()

    for minion_run in minion_runs:

        run_jobs = JobMaster.objects.filter(run=minion_run).filter(running=False)

        for run_job in run_jobs:

            if run_job.job_type.name == "Alignment":

                print("trying to run bwa alignment {} {} {} {}".format(
                    minion_run.id,
                    run_job.id,
                    run_job.reference.id,
                    run_job.last_read
                ))

                run_bwa_alignment.delay(minion_run.id, run_job.id, run_job.reference.id, run_job.last_read)

            if run_job.job_type.name == "Minimap2_trans":
                print ("trying to run transcription alignemnt {} {} {} {}".format(
                    minion_run.id,
                    run_job.id,
                    run_job.reference.id,
                    run_job.last_read
                ))

                run_minimap2_transcriptome.delay(minion_run.id, run_job.id, run_job.reference.id, run_job.last_read)

            if run_job.job_type.name == "Minimap2":

                print("trying to run alignment {} {} {} {}".format(
                    minion_run.id,
                    run_job.id,
                    run_job.reference.id,
                    run_job.last_read
                ))

                run_minimap2_alignment.delay(minion_run.id, run_job.id, run_job.reference.id, run_job.last_read)

            if run_job.job_type.name == "Kraken":

                run_kraken.delay(minion_run.id, run_job.id, run_job.last_read)

            if run_job.job_type.name == "ProcAlign":

                proc_alignment.delay(minion_run.id, run_job.id, run_job.reference.id, run_job.last_read)

            if run_job.job_type.name == "ChanCalc":

                processreads.delay(minion_run.id, run_job.id, run_job.last_read)



@task()
def slow_monitor():

    logger.info('Running slow_monitor celery task.')

    print("Slow Monitor Called")

    testset={}

    cachesiz={}

    minion_runs = MinIONRun.objects.all()

    timediff = utcnow() - timedelta(days=1)

    for minion_run in minion_runs:

        if minion_run.last_entry() >= timediff or minion_run.last_read() >= timediff:

            cachesiz[str(minion_run.id)]=minion_run

    try:
        testset = cache.get('a-unique-key', {})

    except:

        print('a-unique-key not found')

    print('a-unique-key is {}'.format(testset))

    deleted, added = compare_two(testset, cachesiz)

    processrun(deleted,added)

    cache.set('a-unique-key',cachesiz)

    # TODO: We need someway of removing things from the dictionary which aren't still active - otherwise things will
    # persist for ever - so a compare to dictionaries.


def processrun(deleted,added):
    for run in added:
        runinstance = MinIONRun.objects.get(pk=run)
        jobinstance = JobType.objects.get(name="ChanCalc")
        runinstance.active = True
        runinstance.save()

        newjob, created = JobMaster.objects.get_or_create(run=runinstance, job_type=jobinstance)

        if created is True:
            newjob.last_read=0

        newjob.save()

    for run in deleted:
        runinstance = MinIONRun.objects.get(pk=run)
        runinstance.active = False
        runinstance.save()

        jobinstance = JobType.objects.get(name="ChanCalc")
        JobMaster.objects.filter(job_type=jobinstance, run=runinstance).update(complete=True)


def compare_two(newset,cacheset):
    print("deleted keys",newset.keys()-cacheset.keys())

    deleted = newset.keys()-cacheset.keys()

    print("added keys",cacheset.keys() - newset.keys())

    added = cacheset.keys() - newset.keys()

    return deleted,added


@task()
def processreads(runid,id,last_read):
    #print('>>>> Running processreads celery task.')
    print ('running processreads with {} {} {}'.format(runid,id,last_read))
    JobMaster.objects.filter(pk=id).update(running=True)

    fastqs = FastqRead.objects.filter(run_id=runid, id__gt=int(last_read))[:10000]

    chanstore = dict()

    histstore = dict()

    histstore['All reads'] = {}

    barstore = set()

    for fastq in fastqs:
        #print (fastq)
        barstore.add(fastq.barcode)
        #print('>>>> barcode: {}'.format(fastq.barcode))
        if fastq.channel not in chanstore.keys():
            chanstore[fastq.channel]=dict()
            chanstore[fastq.channel]['count']=0
            chanstore[fastq.channel]['length']=0
        chanstore[fastq.channel]['count']+=1
        chanstore[fastq.channel]['length']+=len(fastq.sequence)

        if fastq.barcode.name not in histstore.keys():
            histstore[fastq.barcode.name] = {}

        if fastq.type.name not in histstore[fastq.barcode.name].keys():
            histstore[fastq.barcode.name][fastq.type.name] = {}

        if fastq.type.name not in histstore['All reads'].keys():
            histstore['All reads'][fastq.type.name] = {}

        bin_width = findbin(len(fastq.sequence))

        if bin_width not in histstore[fastq.barcode.name][fastq.type.name].keys():
            histstore[fastq.barcode.name][fastq.type.name][bin_width] = {}

            histstore[fastq.barcode.name][fastq.type.name][bin_width]['count'] = 0
            histstore[fastq.barcode.name][fastq.type.name][bin_width]['length'] = 0

        if bin_width not in histstore['All reads'][fastq.type.name].keys():
            histstore['All reads'][fastq.type.name][bin_width] = {}

            histstore['All reads'][fastq.type.name][bin_width]['count'] = 0
            histstore['All reads'][fastq.type.name][bin_width]['length'] = 0

        histstore[fastq.barcode.name][fastq.type.name][bin_width]['count'] += 1
        histstore[fastq.barcode.name][fastq.type.name][bin_width]['length'] += len(fastq.sequence)

        histstore['All reads'][fastq.type.name][bin_width]['count'] += 1
        histstore['All reads'][fastq.type.name][bin_width]['length'] += len(fastq.sequence)

        ### Adding in the current post save behaviour:

        ipn_obj = fastq

        barcode = ipn_obj.barcode

        barcode_all_reads = ipn_obj.run_id.barcodes.filter(name='All reads').first()

        tm = ipn_obj.start_time

        tm = tm - timedelta(
            minutes=(tm.minute % 1) - 1,
            seconds=tm.second,
            microseconds=tm.microsecond
        )

        obj1, created1 = RunSummaryBarcode.objects.update_or_create(
            run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=barcode_all_reads
        )
        update_sum_stats(obj1, ipn_obj)

        # all reads and barcodes are saved on RunStatisticBarcode
        obj3, created3 = RunStatisticBarcode.objects.update_or_create(
            run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=barcode_all_reads, sample_time=tm
        )
        update_sum_stats(obj3, ipn_obj)

        if barcode is not None and barcode.name != '':
            obj2, created2 = RunSummaryBarcode.objects.update_or_create(
                run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=barcode
            )
            update_sum_stats(obj2, ipn_obj)

            obj3, created3 = RunStatisticBarcode.objects.update_or_create(
                run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=barcode, sample_time=tm
            )
            update_sum_stats(obj3, ipn_obj)

        last_read = fastq.id

    runinstance = MinIONRun.objects.get(pk=runid)

    for chan in chanstore:
        channel, created = ChannelSummary.objects.get_or_create(run_id=runinstance,channel_number=int(chan))
        channel.read_count+=chanstore[chan]['count']
        channel.read_length+=chanstore[chan]['length']
        channel.save()

    for barcodename in histstore.keys():

        barcode = Barcode.objects.filter(run=runinstance, name=barcodename).first()

        for read_type_name in histstore[barcodename].keys():
            read_type = FastqReadType.objects.get(name=read_type_name)
            #print('>>>> {}'.format(read_type_name))

            for hist in histstore[barcodename][read_type_name].keys():
                histogram, created = HistogramSummary.objects.get_or_create(
                    run_id=runinstance,
                    bin_width=int(hist),
                    read_type=read_type,
                    barcode=barcode
                )

                histogram.read_count += histstore[barcodename][read_type_name][hist]["count"]
                histogram.read_length += histstore[barcodename][read_type_name][hist]["length"]
                histogram.save()

    JobMaster.objects.filter(pk=id).update(running=False, last_read=last_read)


def update_sum_stats(obj, ipn_obj):

    if ipn_obj.is_pass:

        obj.pass_length += len(ipn_obj.sequence)

        if len(ipn_obj.sequence) > obj.pass_max_length:
            obj.pass_max_length = len(ipn_obj.sequence)

        if obj.pass_min_length == 0:
            obj.pass_min_length = len(ipn_obj.sequence)

        if len(ipn_obj.sequence) < obj.pass_min_length:
            obj.pass_min_length = len(ipn_obj.sequence)

        obj.pass_count += 1

    obj.total_length += len(ipn_obj.sequence)

    if len(ipn_obj.sequence) > obj.max_length:
        obj.max_length = len(ipn_obj.sequence)

    if obj.min_length == 0:
        obj.min_length = len(ipn_obj.sequence)

    if len(ipn_obj.sequence) < obj.min_length:
        obj.min_length = len(ipn_obj.sequence)

    #
    # channel_presence is a 512 characters length string containing 0
    # for each channel (id between 1 - 512) seen, we set the position on
    # the string to 1. afterwards, we can calculate the number of active
    # channels in a particular minute.
    #
    channel = ipn_obj.channel
    channel_sequence = obj.channel_presence
    channel_sequence_list = list(channel_sequence)
    channel_sequence_list[channel-1] = '1'
    obj.channel_presence = ''.join(channel_sequence_list)

    obj.read_count += 1
    obj.save()


@task()
def proc_alignment(runid,id,reference,last_read):
    JobMaster.objects.filter(pk=id).update(running=True)
    sams = SamStore.objects.filter(run_id=runid,id__gt=int(last_read))[:1000]
    #fp = tempfile.TemporaryFile()
    fp = open('workfile', 'w')
    for sam in sams:
        #print (sam.samline)
        fp.write(sam.samline)
        fp.write('\n')
        last_read=sam.id
    #print (fp.read())
    fp.close()
    #subprocess.run(["samtools", "faidx", "references/hep_ref.fasta"])
    subprocess.run(["samtools", "view", "-bt", "/Volumes/BigElements/human_ref/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa.fai", "workfile", "-o" ,"workfile.bam"])
    #subprocess.run(["rm", "workfile"])
    subprocess.run(["samtools", "sort", "workfile.bam", "-o", "sort_workfile.bam"])
    #subprocess.run(["rm", "workfile.bam"])
    subprocess.run(["samtools", "index", "sort_workfile.bam"])
    stdoutdata = subprocess.getoutput("pysamstats --type variation sort_workfile.bam  --fasta /Volumes/BigElements/human_ref/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa")
    #print("stdoutdata:\r\n " + stdoutdata)
    #subprocess.run(["rm sort_workfile.bam"])

    """
    Then we need to:
        samtools faidx references/hep_ref.fasta
        samtools view -bt references/hep_ref.fasta.fai workfile > workfile.bam
        samtools sort workfile.bam > sort_workfile.bam
        samtools index sort_workfile.bam
        pysamstats --type variation sort_workfile.bam  --fasta references/hep_ref.fasta

    After all that we just (!) parse the pysamstats lines into the existing reference table covering this information.
    """
    JobMaster.objects.filter(pk=id).update(running=False,last_read=last_read)


@task()
def test_task(string, reference):
    REFERENCELOCATION = getattr(settings, "REFERENCELOCATION", None)
    print (reference)
    Reference = ReferenceInfo.objects.get(pk=reference)
    print (type(Reference))
    print (string, REFERENCELOCATION)
    print (Reference)
    print (Reference.minimap2_index_file_location)
    lines = Reference.referencelines.all()
    for line in lines:
        print (line.id)


@task()
def run_minimap2_transcriptome(runid,id,reference,last_read):
    JobMaster.objects.filter(pk=id).update(running=True)
    REFERENCELOCATION = getattr(settings, "REFERENCELOCATION", None)
    Reference = ReferenceInfo.objects.get(pk=reference)
    chromdict = dict()
    chromosomes = Reference.referencelines.all()
    for chromosome in chromosomes:
        chromdict[chromosome.line_name] = chromosome
    minimap2 = Reference.minimap2_index_file_location
    minimap2_ref = os.path.join(REFERENCELOCATION, minimap2)
    # print ("runid:{} last_read:{}".format(runid,last_read))
    fastqs = FastqRead.objects.filter(run_id__id=runid, id__gt=int(last_read))[:1000]
    # print ("fastqs",fastqs)
    read = ''
    fastqdict = dict()
    fastqtypedict = dict()

    # print (len(fastqs))

    for fastq in fastqs:
        read = read + '>{} \r\n{}\r\n'.format(fastq.read_id, fastq.sequence)
        fastqdict[fastq.read_id] = fastq
        fastqtypedict[fastq.read_id] = fastq.type
        last_read = fastq.id
    # print (read)
    cmd = 'minimap2 -x map-ont -t 8 --secondary=no %s -' % (minimap2_ref)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            stdin=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate(input=read.encode("utf-8"))
    status = proc.wait()
    paf = out.decode("utf-8")

    pafdata = paf.splitlines()
    runinstance = MinIONRun.objects.get(pk=runid)

    resultstore = dict()

    for line in pafdata:
        line = line.strip('\n')
        record = line.split('\t')
        # print(record)
        # readid = FastqRead.objects.get(read_id=record[0])
        readid = fastqdict[record[0]]
        typeid = fastqtypedict[record[0]]
        newpaf = PafStore_transcriptome(run=runinstance, read=readid, read_type=typeid)
        newpaf.reference = Reference
        newpaf.qsn = record[0]  # models.CharField(max_length=256)#1	string	Query sequence name
        newpaf.qsl = int(record[1])  # models.IntegerField()#2	int	Query sequence length
        newpaf.qs = int(record[2])  # models.IntegerField()#3	int	Query start (0-based)
        newpaf.qe = int(record[3])  # models.IntegerField()#4	int	Query end (0-based)
        newpaf.rs = record[4]  # models.CharField(max_length=1)#5	char	Relative strand: "+" or "-"
        # newpaf.tsn = record[5] #models.CharField(max_length=256)#6	string	Target sequence name
        newpaf.tsn = chromdict[record[5]]  # models.CharField(max_length=256)#6	string	Target sequence name
        newpaf.tsl = int(record[6])  # models.IntegerField()#7	int	Target sequence length
        newpaf.ts = int(record[7])  # models.IntegerField()#8	int	Target start on original strand (0-based)
        newpaf.te = int(record[8])  # models.IntegerField()#9	int	Target end on original strand (0-based)
        newpaf.nrm = int(record[9])  # models.IntegerField()#10	int	Number of residue matches
        newpaf.abl = int(record[10])  # models.IntegerField()#11	int	Alignment block length
        newpaf.mq = int(record[11])  # models.IntegerField()#12	int	Mapping quality (0-255; 255 for missing)
        newpaf.save()
        if Reference not in resultstore:
            resultstore[Reference] = dict()
        if chromdict[record[5]] not in resultstore[Reference]:
            resultstore[Reference][chromdict[record[5]]] = dict()
        if readid.barcode not in resultstore[Reference][chromdict[record[5]]]:
            resultstore[Reference][chromdict[record[5]]][readid.barcode] = dict()
        if typeid not in resultstore[Reference][chromdict[record[5]]][readid.barcode]:
            resultstore[Reference][chromdict[record[5]]][readid.barcode][typeid] = dict()
        if 'read' not in resultstore[Reference][chromdict[record[5]]][readid.barcode][typeid]:
            resultstore[Reference][chromdict[record[5]]][readid.barcode][typeid]['read'] = set()
            resultstore[Reference][chromdict[record[5]]][readid.barcode][typeid]['length'] = 0
        resultstore[Reference][chromdict[record[5]]][readid.barcode][typeid]['read'].add(record[0])
        resultstore[Reference][chromdict[record[5]]][readid.barcode][typeid]['length'] += int(record[3]) - int(
            record[2]) + 1

    for ref in resultstore:
        for ch in resultstore[ref]:
            for bc in resultstore[ref][ch]:
                for ty in resultstore[ref][ch][bc]:
                    # print (ref,ch,bc,ty,len(resultstore[ref][ch][bc][ty]['read']),resultstore[ref][ch][bc][ty]['length'])
                    summarycov, created2 = PafSummaryCov_transcriptome.objects.update_or_create(
                        run=runinstance,
                        read_type=ty,
                        barcode=bc,
                        reference=ref,
                        chromosome=ch,
                    )
                    summarycov.read_count += len(resultstore[ref][ch][bc][ty]['read'])
                    summarycov.cumu_length += resultstore[ref][ch][bc][ty]['length']
                    summarycov.save()
    # print("!*!*!*!*!*!*!*!*!*!*! ------- running alignment")
    JobMaster.objects.filter(pk=id).update(running=False, last_read=last_read, read_count=F('read_count') + len(fastqs))

@task()
def run_minimap2_alignment(runid, id, reference, last_read):
    JobMaster.objects.filter(pk=id).update(running=True)
    REFERENCELOCATION = getattr(settings, "REFERENCELOCATION", None)
    Reference = ReferenceInfo.objects.get(pk=reference)
    chromdict=dict()
    chromosomes = Reference.referencelines.all()
    for chromosome in chromosomes:
        chromdict[chromosome.line_name]=chromosome
    minimap2 = Reference.minimap2_index_file_location
    minimap2_ref = os.path.join(REFERENCELOCATION, minimap2)
    #print ("runid:{} last_read:{}".format(runid,last_read))
    fastqs = FastqRead.objects.filter(run_id__id=runid, id__gt=int(last_read))[:1000]
    #print ("fastqs",fastqs)
    read = ''
    fastqdict=dict()
    fastqtypedict=dict()

    #print (len(fastqs))

    for fastq in fastqs:
        read = read + '>{} \r\n{}\r\n'.format(fastq.read_id, fastq.sequence)
        fastqdict[fastq.read_id]=fastq
        fastqtypedict[fastq.read_id]=fastq.type
        last_read = fastq.id
    #print (read)
    cmd = 'minimap2 -x map-ont -t 8 --secondary=no %s -' % (minimap2_ref)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            stdin=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate(input=read.encode("utf-8"))
    status = proc.wait()
    paf = out.decode("utf-8")

    pafdata = paf.splitlines()
    runinstance = MinIONRun.objects.get(pk=runid)

    resultstore=dict()

    for line in pafdata:
        line = line.strip('\n')
        record = line.split('\t')
        #print(record)
        #readid = FastqRead.objects.get(read_id=record[0])
        readid=fastqdict[record[0]]
        typeid=fastqtypedict[record[0]]
        newpaf = PafStore(run=runinstance, read=readid, read_type=typeid)
        newpaf.reference = Reference
        newpaf.qsn = record[0] #models.CharField(max_length=256)#1	string	Query sequence name
        newpaf.qsl = int(record[1]) #models.IntegerField()#2	int	Query sequence length
        newpaf.qs  = int(record[2]) #models.IntegerField()#3	int	Query start (0-based)
        newpaf.qe = int(record[3]) #models.IntegerField()#4	int	Query end (0-based)
        newpaf.rs = record[4] #models.CharField(max_length=1)#5	char	Relative strand: "+" or "-"
        #newpaf.tsn = record[5] #models.CharField(max_length=256)#6	string	Target sequence name
        newpaf.tsn = chromdict[record[5]] #models.CharField(max_length=256)#6	string	Target sequence name
        newpaf.tsl = int(record[6]) #models.IntegerField()#7	int	Target sequence length
        newpaf.ts = int(record[7]) #models.IntegerField()#8	int	Target start on original strand (0-based)
        newpaf.te = int(record[8]) #models.IntegerField()#9	int	Target end on original strand (0-based)
        newpaf.nrm = int(record[9]) #models.IntegerField()#10	int	Number of residue matches
        newpaf.abl = int(record[10]) #models.IntegerField()#11	int	Alignment block length
        newpaf.mq = int(record[11]) #models.IntegerField()#12	int	Mapping quality (0-255; 255 for missing)
        newpaf.save()
        if Reference not in resultstore:
            resultstore[Reference]=dict()
        if chromdict[record[5]] not in resultstore[Reference]:
            resultstore[Reference][chromdict[record[5]]]=dict()
        if readid.barcode not in resultstore[Reference][chromdict[record[5]]]:
            resultstore[Reference][chromdict[record[5]]][readid.barcode]=dict()
        if typeid not in resultstore[Reference][chromdict[record[5]]][readid.barcode]:
            resultstore[Reference][chromdict[record[5]]][readid.barcode][typeid]=dict()
        if 'read' not in resultstore[Reference][chromdict[record[5]]][readid.barcode][typeid]:
            resultstore[Reference][chromdict[record[5]]][readid.barcode][typeid]['read']=set()
            resultstore[Reference][chromdict[record[5]]][readid.barcode][typeid]['length']=0
        resultstore[Reference][chromdict[record[5]]][readid.barcode][typeid]['read'].add(record[0])
        resultstore[Reference][chromdict[record[5]]][readid.barcode][typeid]['length']+=int(record[3])-int(record[2])+1

    for ref in resultstore:
        for ch in resultstore[ref]:
            for bc in resultstore[ref][ch]:
                for ty in resultstore[ref][ch][bc]:
                    #print (ref,ch,bc,ty,len(resultstore[ref][ch][bc][ty]['read']),resultstore[ref][ch][bc][ty]['length'])
                    summarycov, created2 = PafSummaryCov.objects.update_or_create(
                        run=runinstance,
                        read_type=ty,
                        barcode=bc,
                        reference=ref,
                        chromosome=ch,
                    )
                    summarycov.read_count += len(resultstore[ref][ch][bc][ty]['read'])
                    summarycov.cumu_length += resultstore[ref][ch][bc][ty]['length']
                    summarycov.save()
    #print("!*!*!*!*!*!*!*!*!*!*! ------- running alignment")
    JobMaster.objects.filter(pk=id).update(running=False, last_read=last_read, read_count=F('read_count')+len(fastqs))


@task()
def run_kraken(runid,id,last_read):
    JobMaster.objects.filter(pk=id).update(running=True)
    fastqs = FastqRead.objects.filter(run_id=runid, id__gt=int(last_read))[:10000]
    krakrun = Kraken()
    # IMPORTANT - NEEDS FIXING FOR BARCODES AND READ TYPES
    if len(fastqs) > 0:
        read=''
        readdict=dict()
        typedict=dict()
        barcodedict=dict()
        for fastq in fastqs:
            read = '{}>{} \r\n{}\r\n'.format(read,fastq,fastq.sequence)
            #print (fastq.id)
            readdict[fastq.read_id]=fastq.id
            last_read=fastq.id
        krakrun.write_seqs(read.encode("utf-8"))
        krakenoutput = krakrun.run().decode('utf-8').split('\n')
        runinstance = MinIONRun.objects.get(pk=runid)
        for line in krakenoutput:
            if len(line) > 0:
                #print (line)
                krakenbits = line.split("\t")


                #print (readdict[krakenbits[1]],krakenbits[1])
                newkrak = MiniKraken(run_id=runid, read_id=int(readdict[krakenbits[1]]))
                newkrak.krakenstatus=str(krakenbits[0])
                newkrak.krakentaxid=int(krakenbits[2])
                newkrak.krakenseqlen=int(krakenbits[3])
                newkrak.krakenlca=str(krakenbits[4])
                newkrak.save()

        #To get the data about the run - specifically which barcodes it has:
        Barcodes = Barcode.objects.filter(run_id=runid)

        #To get the data about the run - specifcally read types:
        Types = FastqReadType.objects.all()

        for Bar in Barcodes:
            for Type in Types:
                #print (Bar,Type)
                #MiniKrak = MiniKraken.objects.filter(run_id=runid).filter(barcode=Bar).filter(read_type=Type)
                MiniKrak = MiniKraken.objects.filter(run=runinstance).filter(read__barcode=Bar).filter(read__type=Type)
                if len(MiniKrak)>0:
                    kraken = ""
                    for krak in MiniKrak:
                        kraken = ("{}{}\t{}\t{}\t{}\t{}\n".format(kraken, krak.krakenstatus, krak.read_id, krak.krakentaxid,
                                                          krak.krakenseqlen, krak.krakenlca))
                    krakrun.write_kraken(kraken.encode("utf-8"))
                    output = krakrun.process().decode('utf-8').split('\n')
                    counter = 0
                    rankdict=dict()
                    rankdict[-1]='Input'
                    for line in output:
                        if len(line)>0:
                            linebits = line.split("\t")
                            #print (linebits)
                            #print (float(linebits[0].lstrip(' ')))
                            parsekrak,created = ParsedKraken.objects.get_or_create(run_id=runid,NCBItaxid=linebits[4],type=Type,barcode=Bar)
                            parsekrak.percentage = float(linebits[0].lstrip(' '))
                            parsekrak.rootreads = int(linebits[1])
                            parsekrak.directreads = int(linebits[2])
                            parsekrak.rank = str(linebits[3])
                            parsekrak.sci_name = str(linebits[5].lstrip(' '))
                            parsekrak.indentation = int((len(linebits[5]) - len(linebits[5].lstrip(' ')))/2)
                            parsekrak.orderin = counter
                            rankdict[int((len(linebits[5]) - len(linebits[5].lstrip(' ')))/2)]=str(linebits[5].lstrip(' '))
                            parsekrak.parent = rankdict[int((len(linebits[5]) - len(linebits[5].lstrip(' ')))/2)-1]
                            parsekrak.save()
                            counter += 1
                    #for krak in MiniKrak:
                        #print (krak)

    krakrun.finish()
    JobMaster.objects.filter(pk=id).update(running=False, last_read=last_read, read_count=F('read_count')+len(fastqs))


class Kraken():
    def __init__(self):
        self.tmpfile = tempfile.NamedTemporaryFile(suffix=".fa")
        self.krakenfile = tempfile.NamedTemporaryFile(suffix=".out")
        #db

    def write_seqs(self, seqs):
        #self.tmpfile.write("\n".join(seqs))
        self.tmpfile.write(seqs)


    def write_kraken(self,kraken):
        #print "printing the kraken"
        #print "".join(kraken)
        self.krakenfile.write(kraken)

    def read_kraken(self):
        for line in self.krakenfile:
            print (line)

    def process(self):
        print ('kraken-report --db /Volumes/SSD/kraken/minikraken_20141208 '+self.krakenfile.name)
        p1 = subprocess.Popen('kraken-report --db /Volumes/SSD/kraken/minikraken_20141208 '+self.krakenfile.name, shell=True, stdout=subprocess.PIPE)
        (out,err) =  p1.communicate()
        return out

    def process2(self):
        print ('kraken-mpa-report --db /Volumes/SSD/kraken/minikraken_20141208 '+self.krakenfile.name)
        p1 = subprocess.Popen('kraken-mpa-report --db /Volumes/SSD/kraken/minikraken_20141208 '+self.krakenfile.name, shell=True, stdout=subprocess.PIPE)
        (out,err) =  p1.communicate()
        return out

    def run(self):
        print ('kraken --quick --db /Volumes/SSD/kraken/minikraken_20141208 --fasta-input --preload  '+self.tmpfile.name) #+'  #| kraken-translate --db /Volumes/SSD/kraken/minikraken_20141208/ $_'
        #p1 = subprocess.Popen('kraken --quick --db /Volumes/SSD/kraken/minikraken_20141208 --fasta-input --preload  '+self.tmpfile.name+'  | kraken-translate --db /Volumes/SSD/kraken/minikraken_20141208/ $_', shell=True, stdout=subprocess.PIPE)
        p1 = subprocess.Popen('kraken --quick --db /Volumes/SSD/kraken/minikraken_20141208 --fasta-input --preload  '+self.tmpfile.name+'  ', shell=True, stdout=subprocess.PIPE)
        #print 'minimap -Sw5 -L100 -t8 '+self.tmpfile.name+' '+self.tmpfile.name+' | miniasm -f '+self.tmpfile.name+' - | awk \'/^S/{print \">\"$2\"\\n\"$3}\' | fold '
        #p1 = subprocess.Popen('minimap -Sw5 -L100 -t8 '+self.tmpfile.name+' '+self.tmpfile.name+' | miniasm -f '+self.tmpfile.name+' - | awk \'/^S/{print \">\"$2\"\\n\"$3}\' | fold ', shell=True, stdout=subprocess.PIPE)
        (out, err) = p1.communicate()
        return out

    def finish(self):
        self.tmpfile.close()
        self.krakenfile.close()


@task()
def run_bwa_alignment(runid, job_id, referenceinfo_id, last_read):

    print('---> Inside bwa alignment')

    job = JobMaster.objects.get(pk=job_id)

    if job.running is True:
        print('Job is already running. No further action required.')
        return

    job.running = True
    job.save()

    print('---> Job was not running')

    # JobMaster.objects.filter(pk=job_id).update(running=True)

    #print (runid,id,reference,last_read)
    fastq_list = FastqRead.objects.filter(
        run_id=runid,
        id__gt=last_read
    )[:1000]

    referenceinfo = ReferenceInfo.objects.get(pk=referenceinfo_id)

    for fastq in fastq_list:
        #print (fastq.sequence)
        bwaindex = referenceinfo.filename
        read = '>{} \r\n{}\r\n'.format(fastq, fastq.sequence)
        cmd = 'bwa mem -x ont2d %s -' % (bwaindex)
        #print (cmd)
        #print (read)
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdin=subprocess.PIPE,
            shell=True
        )

        (out, err) = proc.communicate(input=read.encode("utf-8"))
        status = proc.wait()
        sam = out.decode("utf-8")

        print('---> sam')

        samdata = sam.splitlines()

        print(samdata)

        for line in samdata:
            if not line.startswith('@'):
                line = line.strip('\n')
                record = line.split('\t')
                if record[2] != '*':
                    runinstance = MinIONRun.objects.get(pk=runid)
                    newsam = SamStore(run_id=runinstance, read_id=fastq)
                    newsam.samline = line
                    newsam.save()
                    last_read=fastq.id

    job.running = False
    job.last_read = last_read
    job.save()


@task
def updateReadNamesOnRedis():
    print ('>>> running updateReadNamesOnRedis')
    r = redis.StrictRedis(host='localhost', port=6379, db=0)

    runs = MinIONRun.objects.all()

    for run in runs:
        print ('>>> run: {}'.format(run.id))

        reads = FastqRead.objects.filter(run_id=run).order_by('id')

        paginator = Paginator(reads, 25)

        key = 'run.{}.reads.number_pages'.format(run.id)

        r.set(key, paginator.num_pages)

        for page in range(1, paginator.num_pages + 1):
            print ('>>> run {}, page {} of {}'.format(run.id, page, paginator.num_pages))

            key = 'run.{}.reads.page.{}'.format(run.id, page)

            result = paginator.page(page)

            result2 = set()

            for key in result:
                result2.add(key.read_id)

            r.set(key, json.dumps(list(result2)))


@task
def delete_runs():

    MinIONRun.objects.filter(to_delete=True).delete()


@task
def send_messages():

    new_messages = Message.objects.filter(delivered_date=None)

    for new_message in new_messages:

        #print('Sending message: {}'.format(new_message))

        message_sent = False

        if new_message.recipient.extendedopts.email:

            try:

                if new_message.sender:

                    sender = new_message.sender.email

                else:

                    sender = 'contact@minotour.org'

                send_mail(
                    new_message.title,
                    new_message.content,
                    sender,
                    [new_message.recipient.email],
                    fail_silently=False,
                )

                message_sent = True

            except MailgunAPIError as e:
                print(e)

        if new_message.recipient.extendedopts.tweet \
            and new_message.recipient.extendedopts.twitterhandle != '':

            TWITTOKEN = settings.TWITTOKEN
            TWITTOKEN_SECRET = settings.TWITTOKEN_SECRET
            TWITCONSUMER_KEY = settings.TWITCONSUMER_KEY
            TWITCONSUMER_SECRET = settings.TWITCONSUMER_SECRET

            t = Twitter(
                auth=OAuth(TWITTOKEN, TWITTOKEN_SECRET, TWITCONSUMER_KEY, TWITCONSUMER_SECRET)
            )

            t.direct_messages.new(
                user=new_message.recipient.extendedopts.twitterhandle,
                text=new_message.title
            )
            #status = '@{} {}'.format(new_message.recipient.extendedopts.twitterhandle,new_message.title)
            #t.statuses.update(
            #    status=status
            #)

            message_sent = True

        if message_sent:
            print('inside message_sent')
            new_message.delivered_date = utcnow()
            new_message.save()
