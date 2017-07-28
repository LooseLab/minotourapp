from __future__ import absolute_import, unicode_literals
from celery import task
from reads.models import MinIONRun
from reads.models import JobMaster
from reads.models import FastqRead
from reads.models import SamStore
from reads.models import ChannelSummary
from reads.models import Job
from datetime import datetime, timedelta
from django.db.models import Q
import subprocess
import tempfile
from django.core.cache import cache

@task()
def run_monitor():
    # Do something...
    print ("Rapid Monitor Called")
    minion_runs = MinIONRun.objects.filter(Q(reads__created_date__gte=datetime.now() - timedelta(days=1)) | Q(
            RunStats__created_date__gte=datetime.now() - timedelta(days=1))).distinct()
    print(minion_runs)
    print(len(minion_runs))
    for minion_run in minion_runs:
        print ("checking jobs for {}".format(minion_run))
        print (minion_run.run_id)
        run_jobs = JobMaster.objects.filter(run_id=minion_run.id)
        for run_job in run_jobs:
            print (type(run_job.job_name))
            if str(run_job.job_name)=="Alignment" and run_job.running is False:
                print ("trying to run alignment")
                run_alignment.delay(minion_run.id,run_job.id,run_job.var1,run_job.var2)
            if run_job.running is True:
                print ("{} is already running!".format(run_job.job_name) )

@task()
def slow_monitor():
    # Do something...
    print ("Slow Monitor Called")
    testset={}
    cachesiz={}
    minion_runs = MinIONRun.objects.filter(Q(reads__created_date__gte=datetime.now() - timedelta(days=1)) | Q(
            RunStats__created_date__gte=datetime.now() - timedelta(days=1))).distinct()
    print(minion_runs)
    print(len(minion_runs))
    for minion_run in minion_runs:
        print ("checking jobs for {}".format(minion_run))
        print (minion_run.run_id)
        cachesiz[str(minion_run.run_id)]=minion_run
        run_jobs = JobMaster.objects.filter(run_id=minion_run.id)
        for run_job in run_jobs:
            print (type(run_job.job_name))
            if str(run_job.job_name)=="ProcAlign" and run_job.running is False:
                print ("trying to process alignment")
                proc_alignment.delay(minion_run.id, run_job.id, run_job.var1, run_job.var2)
            if str(run_job.job_name)=="ChanCalc" and run_job.running is False:
                print ("ChannelCalc")
                calculate_channels.delay(minion_run.id, run_job.id, run_job.var1, run_job.var2)
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
        newjob = JobMaster(run_id=runinstance,job_name=jobinstance,var2=0)
        newjob.save()
    for run in deleted:
        runinstance = MinIONRun.objects.get(run_id=run)
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
def calculate_channels(runid,id,var1,last_read):
    JobMaster.objects.filter(pk=id).update(running=True)
    fastqs = FastqRead.objects.filter(run_id=runid, id__gt=int(last_read))[:10000]
    tempstore=dict()
    for fastq in fastqs:
        #print (fastq)
        if fastq.channel not in tempstore.keys():
            tempstore[fastq.channel]=dict()
            tempstore[fastq.channel]['count']=0
            tempstore[fastq.channel]['length']=0
        tempstore[fastq.channel]['count']+=1
        tempstore[fastq.channel]['length']+=len(fastq.sequence)
        last_read = fastq.id
    #print (tempstore)
    runinstance = MinIONRun.objects.get(pk=runid)
    for chan in tempstore:
        #print (chan)
        channel, created = ChannelSummary.objects.get_or_create(run_id=runinstance,channel_number=int(chan),read_count=0,read_length=0)
        #print (tempstore[chan]['count'])
        channel.read_count+=tempstore[chan]['count']
        channel.read_length+=tempstore[chan]['length']
        channel.save()
    JobMaster.objects.filter(pk=id).update(running=False, var2=last_read)


@task()
def proc_alignment(runid,id,reference,last_read):
    JobMaster.objects.filter(pk=id).update(running=True)
    sams = SamStore.objects.filter(run_id=runid,id__gt=int(last_read))[:250]
    #fp = tempfile.TemporaryFile()
    fp = open('workfile', 'w')
    for sam in sams:
        #print (sam.samline)
        fp.write(sam.samline)
        fp.write('\n')
        last_read=sam.id
    #print (fp.read())
    fp.close()
    subprocess.run(["samtools", "faidx", "references/hep_ref.fasta"])
    subprocess.run(["samtools", "view", "-bt", "references/hep_ref.fasta.fai", "workfile", "-o" ,"workfile.bam"])
    #subprocess.run(["rm", "workfile"])
    subprocess.run(["samtools", "sort", "workfile.bam", "-o", "sort_workfile.bam"])
    #subprocess.run(["rm", "workfile.bam"])
    subprocess.run(["samtools", "index", "sort_workfile.bam"])
    stdoutdata = subprocess.getoutput("pysamstats --type variation sort_workfile.bam  --fasta references/hep_ref.fasta")
    print("stdoutdata:\r\n " + stdoutdata)
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
def run_alignment(runid,id,reference,last_read):
    JobMaster.objects.filter(pk=id).update(running=True)
    print (runid,id,reference,last_read)
    fastqs = FastqRead.objects.filter(run_id=runid,id__gt=int(last_read))[:250]
    for fastq in fastqs:
        #print (fastq.sequence)
        bwaindex = 'references/hep_ref.fasta'
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

