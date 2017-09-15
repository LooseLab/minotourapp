from __future__ import absolute_import, unicode_literals
from celery import task
from celery.utils.log import get_task_logger

from reads.models import MinIONRun, FastqReadType
from reads.models import JobMaster
from reads.models import FastqRead
from alignment.models import SamStore
from reads.models import ChannelSummary
from reads.models import HistogramSummary
from reads.models import Job
from reads.models import UserOptions
from reads.models import Barcode
from reference.models import ReferenceInfo
from reference.models import ReferenceLine
from alignment.models import PafStore
from datetime import datetime, timedelta
from django.db.models import Q
import subprocess
import tempfile
from django.core.cache import cache
from django.contrib.auth.models import User
from django.conf import settings
import pytz
from twitter import *
import os

logger = get_task_logger(__name__)


def send_tweet(message):
    TWITTOKEN = getattr(settings, "TWITTOKEN", None)
    TWITTOKEN_SECRET = getattr(settings, "TWITTOKEN_SECRET", None)
    TWITCONSUMER_KEY = getattr(settings, "TWITCONSUMER_KEY", None)
    TWITCONSUMER_SECRET = getattr(settings, "TWITCONSUMER_SECRET", None)
    t = Twitter(
        auth=OAuth(TWITTOKEN, TWITTOKEN_SECRET, TWITCONSUMER_KEY, TWITCONSUMER_SECRET))
    print ("Gonna try and send a tweet")
    t.statuses.update(
        status=message)

def utcnow():
    return datetime.now(tz=pytz.utc)


@task()
def run_monitor():
    logger.info('Running run_monitor celery task.')
    # Do something...
    print ("Rapid Monitor Called")
    #minion_runs = MinIONRun.objects.filter(Q(reads__created_date__gte=datetime.utcnow() - timedelta(days=3)) | Q(
    #        RunStats__created_date__gte=datetime.now() - timedelta(days=3))).distinct()
    minion_runs = MinIONRun.objects.filter(active=True).distinct()
    #print(minion_runs)
    #print(len(minion_runs))
    for minion_run in minion_runs:
        #print ("Run Monitor jobs for {}".format(minion_run))
        #print (minion_run.run_id)
        #print (minion_run.owner)
        #SendUserMessage.delay(minion_run.id,"a test","testing")
        run_jobs = JobMaster.objects.filter(run_id=minion_run.id)
        for run_job in run_jobs:
            #print (type(run_job.job_name))
            if str(run_job.job_name)=="Alignment" and run_job.running is False:
                #print ("trying to run alignment")
                run_alignment.delay(minion_run.id,run_job.id,run_job.var1,run_job.var2)
            if str(run_job.job_name)=="Minimap2" and run_job.running is False:
                #print ("trying to run alignment")
                run_minimap2.delay(minion_run.id,run_job.id,run_job.var1.id,run_job.var2)
            if run_job.running is True:
                #print ("{} is already running!".format(run_job.job_name) )
                pass

@task()
def slow_monitor():
    logger.info('Running slow_monitor celery task.')
    # Do something...
    print ("Slow Monitor Called")
    testset={}
    cachesiz={}
    minion_runs = MinIONRun.objects.filter(Q(reads__created_date__gte=utcnow() - timedelta(days=10)) | Q(
            RunStats__created_date__gte=utcnow() - timedelta(days=10))).distinct()
    #minion_runs = MinIONRun.objects.all()
    print(minion_runs)
    #print(len(minion_runs))
    for minion_run in minion_runs:
        #print ("checking jobs for {}".format(minion_run))
        #print (minion_run.run_id)
        cachesiz[str(minion_run.run_id)]=minion_run
        run_jobs = JobMaster.objects.filter(run_id=minion_run.id)
        for run_job in run_jobs:
            #print('>>> run_jobs')
            #print(run_job.job_name)
            #print(run_job.running)

            if str(run_job.job_name)=="ProcAlign" and run_job.running is False:
                #print ("trying to process alignment")
                proc_alignment.delay(minion_run.id, run_job.id, run_job.var1, run_job.var2)
            if str(run_job.job_name)=="ChanCalc" and run_job.running is False:
                #print ("ChannelCalc")
                processreads.delay(minion_run.id, run_job.id, run_job.var1, run_job.var2)
    try:
        testset = cache.get('a-unique-key', {})
    except:
        print('a-unique-key not found')
    print('a-unique-key is {}'.format(testset))
    deleted,added = compare_two(testset,cachesiz)
    processrun(deleted,added)
    cache.set('a-unique-key',cachesiz)
    ### We need someway of removing things from the dictionary which aren't still active - otherwise things will persist for ever - so a compare to dictionaries.

def processrun(deleted,added):
    for run in added:
        runinstance = MinIONRun.objects.get(run_id=run)
        jobinstance = Job.objects.get(jobname="ChanCalc")
        runinstance.active=True
        runinstance.save()
        #newjob = JobMaster(run_id=runinstance,job_name=jobinstance,var2=0)
        newjob,created = JobMaster.objects.get_or_create(run_id=runinstance, job_name=jobinstance)
        if created is True:
            newjob.var2=0
        newjob.save()
    for run in deleted:
        runinstance = MinIONRun.objects.get(run_id=run)
        runinstance.active=False
        runinstance.save()
        jobinstance = Job.objects.get(jobname="ChanCalc")
        jobrecord=JobMaster.objects.filter(job_name=jobinstance,run_id=runinstance).update(complete=True)

def compare_two(newset,cacheset):
    # Do somehing ....
    #print ("newset",newset.keys())
    #print ("cached",cacheset.keys())
    print ("deleted keys",newset.keys()-cacheset.keys())
    deleted = newset.keys()-cacheset.keys()
    print("added keys",cacheset.keys() - newset.keys())
    added = cacheset.keys() - newset.keys()
    return deleted,added

@task()
def SendUserMessage(runid,messagetype,messagestring):
    print ("Message Sending Initiated")
    print ("looking for {}".format(runid))
    runinstance = MinIONRun.objects.get(id=runid)
    print ("and now for {}".format(runinstance.owner))
    UserObject = User.objects.get(username=runinstance.owner)
    print (UserObject.extendedopts.tweet)
    print (runinstance,UserObject)


@task()
def processreads(runid,id,var1,last_read):
    #print('>>>> Running processreads celery task.')

    JobMaster.objects.filter(pk=id).update(running=True)
    fastqs = FastqRead.objects.filter(run_id=runid, id__gt=int(last_read))[:10000]
    chanstore=dict()

    histstore=dict()
    histstore['All reads'] = {}

    barstore=set()

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

        last_read = fastq.id

    #print (tempstore)
    runinstance = MinIONRun.objects.get(pk=runid)
    #for barcode in barstore:
    #    result, created = Barcode.objects.get_or_create(run=runinstance, name=barcode)

    for chan in chanstore:
        #print (chan)
        channel, created = ChannelSummary.objects.get_or_create(run_id=runinstance,channel_number=int(chan))
        #print ('created',created)
        #print (channel)
        #print (tempstore[chan]['count'])
        channel.read_count+=chanstore[chan]['count']
        channel.read_length+=chanstore[chan]['length']
        channel.save()

    #print('>>>> {}'.format(histstore))

    for barcodename in histstore.keys():
        #print('>>>> {}'.format(barcodename))

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

    JobMaster.objects.filter(pk=id).update(running=False, var2=last_read)


def findbin(x,bin_width=900):
    return ((x - x % bin_width) / bin_width)

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
    JobMaster.objects.filter(pk=id).update(running=False,var2=last_read)

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
def run_minimap2(runid,id,reference,last_read):
    JobMaster.objects.filter(pk=id).update(running=True)
    REFERENCELOCATION = getattr(settings, "REFERENCELOCATION", None)
    Reference = ReferenceInfo.objects.get(pk=reference)
    chromdict=dict()
    chromosomes = Reference.referencelines.all()
    for chromosome in chromosomes:
        chromdict[chromosome.line_name]=chromosome
    minimap2 = Reference.minimap2_index_file_location
    minimap2_ref = os.path.join(REFERENCELOCATION, minimap2)
    fastqs = FastqRead.objects.filter(run_id=runid, id__gt=int(last_read))[:250]
    read = ''
    fastqdict=dict()

    for fastq in fastqs:
        read = read + '>{} \r\n{}\r\n'.format(fastq.read_id, fastq.sequence)
        fastqdict[fastq.read_id]=fastq
        last_read = fastq.id
    cmd = 'minimap2 -x map-ont -t 8 -N 0 %s -' % (minimap2_ref)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            stdin=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate(input=read.encode("utf-8"))
    status = proc.wait()
    paf = out.decode("utf-8")
    pafdata = paf.splitlines()
    runinstance = MinIONRun.objects.get(pk=runid)

    for line in pafdata:
        line = line.strip('\n')
        record = line.split('\t')
        #print(record)
        #readid = FastqRead.objects.get(read_id=record[0])
        readid=fastqdict[record[0]]
        newpaf = PafStore(run=runinstance, read=readid)
        newpaf.reference_id = ReferenceInfo.objects.get(pk=reference)
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

    JobMaster.objects.filter(pk=id).update(running=False, var2=last_read)
    

@task()
def run_alignment(runid,id,reference,last_read):
    JobMaster.objects.filter(pk=id).update(running=True)
    #print (runid,id,reference,last_read)
    fastqs = FastqRead.objects.filter(run_id=runid,id__gt=int(last_read))[:1000]
    for fastq in fastqs:
        #print (fastq.sequence)
        bwaindex = 'Human'
        read = '>{} \r\n{}\r\n'.format(fastq,fastq.sequence)
        cmd = 'bwa mem -x ont2d %s -' % (bwaindex)
        #print (cmd)
        #print (read)
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(input=read.encode("utf-8"))
        status = proc.wait()
        sam = out.decode("utf-8")
        #print (sam)
        samdata = sam.splitlines()
        for line in samdata:
            if not line.startswith('@'):
                line = line.strip('\n')
                record = line.split('\t')
                if record[2] != '*':
                    runinstance = MinIONRun.objects.get(pk=runid)
                    #print (runinstance)
                    newsam = SamStore(run_id=runinstance, read_id=fastq)
                    newsam.samline = line
                    #newsam.qname = record[0]
                    #print (fastq.quality)
                    #newsam.qualscores = fastq.quality ## Need to check if this is correct - not convinced!
                    #newsam.flag = int(record[1])
                    #newsam.rname = record[2]
                    #newsam.refid = 1 ## Needs to be fixed going forwards
                    #newsam.pos = int(record[3])
                    #newsam.mapq = int(record[4])
                    #newsam.cigar = record[5]
                    #newsam.rnext = record[6]
                    #newsam.pnext = record[7]
                    #newsam.tlen = int(record[8])
                    #newsam.seq = record[9]
                    #newsam.qual = record[10]
                    #newsam.nm = record[11]
                    #newsam.md = record[12]
                    #newsam.ass = record[13]
                    #newsam.xs = record[14]
                    newsam.save()
                    last_read=fastq.id
                    #print (last_read)
    JobMaster.objects.filter(pk=id).update(running=False,var2=last_read)
