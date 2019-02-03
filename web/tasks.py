from __future__ import absolute_import, unicode_literals

import json
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from datetime import datetime, timedelta

import numpy as np
import pytz
import redis
from celery import task
from celery.utils.log import get_task_logger
from django.conf import settings
from django.core.mail import send_mail
from django.core.paginator import Paginator
from django.db.models import Max,Min
from django_mailgun import MailgunAPIError
from twitter import *

from assembly.models import GfaStore, GfaSummary
from centrifuge import centrifuge
from communication.utils import *
from jobs.models import JobMaster
from reads.models import Barcode, FastqRead, Run, FlowcellSummaryBarcode, Flowcell, MinIONRunStatus
from web.tasks_chancalc import chancalc
from .tasks_alignment import run_minimap2_alignment
from centrifuge.sankey import calculate_sankey
from centrifuge.tasks import output_parser

import pandas as pd
import gzip


logger = get_task_logger(__name__)


def utcnow():
    return datetime.now(tz=pytz.utc)


def findbin(x, bin_width=900):
    return (x - x % bin_width) / bin_width


def getn50(lens):
    h = sum(lens)/2
    t = 0
    for l in lens:
        t += l
        if t >= h:
            return l


@task()
def run_monitor():

    logger.info('--------------------------------')
    logger.info('Running run_monitor celery task.')
    logger.info('--------------------------------')

    flowcell_list = Flowcell.objects.filter(is_active=True)

    for flowcell in flowcell_list:

        flowcell_job_list = JobMaster.objects.filter(flowcell=flowcell).filter(running=False, complete=False)

        for flowcell_job in flowcell_job_list:

            flowcell_job.running = True
            flowcell_job.save()

            if flowcell_job.job_type.name == "Minimap2":

                logger.info("Sending task minimap2 to server - Flowcell id: {}, job_master id: {}".format(
                    flowcell.id,
                    flowcell_job.id,
                ))

                run_minimap2_alignment.delay(
                    flowcell.id,
                    flowcell_job.id,
                    flowcell_job.reference.id,
                    flowcell_job.last_read
                )

            if flowcell_job.job_type.name == "ChanCalc":

                logger.info("Sending task chancalc to server - Flowcell id: {}, job_master id: {}".format(
                    flowcell.id,
                    flowcell_job.id,
                ))

                chancalc.delay(flowcell.id, flowcell_job.id, flowcell_job.last_read)

            if flowcell_job.job_type.name == "Assembly":

                logger.info("Sending task assembly to server - Flowcell id: {}, job_master id: {}".format(
                    flowcell.id,
                    flowcell_job.id,
                ))

                inputtype = "flowcell"

                run_minimap2_assembly.delay(flowcell_job.id)

            if flowcell_job.job_type.name == "Metagenomics":

                logger.info("Sending task metagenomics to server - Flowcell id: {}, job_master id: {}".format(
                    flowcell.id,
                    flowcell_job.id,
                ))

                run_centrifuge.delay(flowcell_job.id)

            if flowcell_job.job_type.name == "UpdateFlowcellDetails":

                logger.info("Sending task updateflowcelldetails to server - Flowcell id: {}, job_master id: {}".format(
                    flowcell.id,
                    flowcell_job.id,
                ))

                update_flowcell_details.delay(flowcell_job.id)
            if flowcell_job.job_type.name == "CalculateSankey":
                logger.info("Sending task CalculateSankey - Flowcell id: {}, job_master id: {}".format(
                    flowcell.id,
                    flowcell_job.id,
                ))
                run_sankey(flowcell_job.id)

            if flowcell_job.job_type.name == "Parser":
                logger.info("Sending task Parser - Flowcell id: {}, job_master id: {}".format(
                    flowcell.id,
                    flowcell_job.id,
                ))
                output_parser.delay(flowcell_job.id)

@task()
def run_sankey(flowcell_job_id):
    """
    Calculate sankeys for a flowcell upon demand of the user
    :param flowcell_job_id: The pk id of this task
    :author: Rory
    :return:
    """
    job_master = JobMaster.objects.get(pk=flowcell_job_id)
    logger.info("Flowcell id: {} - Starting Sankey calculation task".format(job_master.flowcell.id))
    calculate_sankey(flowcell_job_id)


@task()
def run_centrifuge(flowcell_job_id):

    job_master = JobMaster.objects.get(pk=flowcell_job_id)

    logger.info("Flowcell id: {} - Starting centrifuge task".format(job_master.flowcell.id))

    centrifuge.run_centrifuge(flowcell_job_id)


@task()
def export_reads(runid,id,tmp,last_read,inputtype):
    """
    Function to write out fastq data to known temporary file, compress it,
    and enable it for export. File will be set to self destruct after 7 days or
    1 day after downloading.
    Function should check to see if the file already exists somewhere.
    Also should split data out by barcode and - preferably read type.
    :param runid:
    :param id:
    :param tmp:
    :param last_read:
    :param inputtype:
    :return:
    """
    JobMaster.objects.filter(pk=id).update(running=True)
    runidset = get_runidset(runid,inputtype)
    barcodeset,barcoded=get_barcodes(runidset)
    # fetch all the reads from the run
    fastqs = FastqRead.objects.filter(run_id__id__in=runidset)
    dirpath = tempfile.mkdtemp()
    filedict={}
    if barcoded==True:
        for barcode in barcodeset:
            filehandle = os.path.join(dirpath,barcode+"_pass"+".fastq")
            filedict.setdefault(filehandle, open(filehandle, 'w'))
            filehandle = os.path.join(dirpath,barcode+"_fail"+".fastq")
            filedict.setdefault(filehandle, open(filehandle, 'w'))
            #files_dict[filename].write("test")
    else:
        filehandle = os.path.join(dirpath,"AllReads_pass.fastq")
        filedict.setdefault(filehandle, open(filehandle, 'w'))
        filehandle = os.path.join(dirpath,"AllReads_fail.fastq")
        filedict.setdefault(filehandle, open(filehandle, 'w'))

    if barcoded==True:
        for fastq in fastqs:
            if fastq.is_pass:
                filehandle = os.path.join(dirpath,fastq.barcode.name+"_pass"+".fastq")
            else:
                filehandle = os.path.join(dirpath, fastq.barcode.name + "_fail" + ".fastq")
            filedict[filehandle].write(format_read(fastq))
    else:
        for fastq in fastqs:
            if fastq.is_pass:
                filehandle = os.path.join(dirpath,"AllReads_pass"+".fastq")
            else:
                filehandle = os.path.join(dirpath, "AllReads_fail" + ".fastq")
            filedict[filehandle].write(format_read(fastq))
    archive_zip = zipfile.ZipFile(os.path.join(dirpath,'archive.zip'), 'w')
    for folder, subfolders, files in os.walk(dirpath):
        for file in files:
            if file.endswith('.fastq'):
                archive_zip.write(os.path.join(folder, file),
                                  os.path.relpath(os.path.join(folder, file), dirpath),
                                  compress_type=zipfile.ZIP_DEFLATED)
                try:
                    os.remove(os.path.join(folder, file))
                except OSError:
                    pass
    archive_zip.close()
    JobMaster.objects.filter(pk=id).update(tempfile_name=dirpath,complete=1)
    later = datetime.utcnow() + timedelta(minutes=5)
    delete_folder.apply_async((id), eta=later)
    send_message([JobMaster.objects.filter(pk=id).flowcell.owner], "Read Export Complete",
                 "Your reads have been exported and can be found here: {}.".format(dirpath))
    JobMaster.objects.filter(pk=id).update(running=False)

@task
def delete_folder(id):
    """
    This task will delete a specific folder based on the tempfile_name in a specific job.
    It will then delete the jobmaster.
    :param id:
    :return:
    """
    jobtoprocess = JobMaster.objects.get(pk=id)
    dirpath = jobtoprocess.tempfile_name
    shutil.rmtree(dirpath)
    jobtoprocess.delete()


def get_barcodes(runidset):
    """
    Given a set of runs, returns all the barcodes within them.
    :param runidset:
    :return:
    """
    barcodeset = set()
    barcodes = Barcode.objects.filter(run_id__in=runidset)
    for barcode in barcodes:
        if barcode.name != "All reads" and barcode.name != "No barcode":
            barcodeset.add(barcode.name)
    barcoded=False
    if len(barcodeset)>0:
        barcoded=True
    return (barcodeset,barcoded)


def format_read(fastq):
    """
    This function takes a fastq database object and returns a formatted string for writing
    to a file
    :param fastq: fastq object from the database
    :return: string formatted version of fastq sequence to write to a file.
    """
    lineheader = ">"+str(fastq)
    if len(fastq.run_id.run_id) > 0:
        lineheader = lineheader + " runid={}".format(fastq.run_id.run_id)
    if len(str(fastq.read)) > 0:
        lineheader = lineheader + " read={}".format(fastq.read)
    if len(str(fastq.channel)) > 0:
        lineheader = lineheader + " channel={}".format(fastq.channel)
    ##This output format needs checking
    if len(str(fastq.start_time)) > 0:
        lineheader = lineheader + " start_time={}".format(fastq.start_time.replace(tzinfo=pytz.UTC).isoformat())
    if fastq.barcode.name != "No barcode":
        lineheader = lineheader + " barcode={}".format(fastq.barcode.name)
    return ("{}\n{}\n+\n{}\n".format(lineheader,fastq.sequence,fastq.quality))


@task()
def clean_up_assembly_files(runid,id,tmp):
    """
    This task will automatically delete a temporary file at some given time period after creation.
    :param runid:
    :param id:
    :param tmp:
    :return:
    """
    os.remove(tmp)
    JobMaster.objects.filter(pk=id).update(tempfile_name="")


def get_runidset(runid, inputtype):

    run_set = set()

    if inputtype == "flowcell":

        flowcell = Flowcell.objects.get(pk=runid)

        for run in flowcell.runs.all():

            run_set.add(run.id)

    else:

        run_set.add(runid)

    return run_set

def callfetchreads_assem(runs,chunk_size,last_read):
    fastasm=list()
    #fastqs_list = list()
    fastq_df_barcode = pd.DataFrame()
    while True:
        #reads, last_read, read_count, fastasmchunk, fastqs_list_chunk = fetchreads_cent(runs, chunk_size, last_read)
        reads, last_read, read_count, fastasmchunk = fetchreads_assem(runs, chunk_size, last_read)
        fastasm+=fastasmchunk
        #fastqs_list += fastqs_list_chunk
        fastq_df_barcode = fastq_df_barcode.append(reads)
        if len(fastq_df_barcode)>=chunk_size or len(reads)==0:
            break
    read_count = len(fastq_df_barcode)
    return fastq_df_barcode.reset_index(),last_read,read_count,fastasm #,fastqs_list


def fetchreads_assem(runs,chunk_size,last_read):
    countsdict=dict()
    fastasm = list()
    #fastqs_list = list()
    for run in runs:
        #fastqs = FastqRead.objects.filter(run=run, id__gt=int(last_read)).first()
        fastqs = FastqRead.objects.values_list('id').filter(run=run, id__gt=int(last_read)).first()
        if fastqs:
            countsdict[fastqs[0]] = run
    count = 1
    fastq_df_barcode = pd.DataFrame()
    if len(countsdict)>1:
        for entry in sorted(countsdict):
            if count < len(countsdict):
                fastqs = FastqRead.objects.filter(run=countsdict[entry],
                                              id__gt=int(last_read),
                                                id__lt=int(sorted(countsdict)[count]),
                                                         )[:chunk_size]
                fastq_df_barcode = fastq_df_barcode.append(pd.DataFrame.from_records(
                    fastqs.values("read_id", "barcode_name", "sequence", "id")))

                # Create a list of the fastQRead objects we will use in this objects
                fastasm += list(fastqs)
                # Create a list of 8000 tuples each with these 4 objects in it
                #fastqs_list += fastqs.values_list('read_id', 'barcode_name', 'sequence', 'id')
                last_read = fastq_df_barcode.id.max()
                if len(fastq_df_barcode) >= chunk_size:
                    break
            count += 1
    elif len(countsdict)==1:
        """This is risky and it breaks the logic - we end up grabbing reads"""
        mykey = list(countsdict.keys())[0]
        fastqs = FastqRead.objects.filter(run=countsdict[mykey],
                                              id__gt=int(last_read),)[:chunk_size]
        fastq_df_barcode = fastq_df_barcode.append(pd.DataFrame.from_records(
            fastqs.values("read_id", "barcode_name", "sequence", "id")))
        # Create a list of the fastQRead objects we will use in this objects
        fastasm += list(fastqs)
        # Create a list of 8000 tuples each with these 4 objects in it
        #fastqs_list += fastqs.values_list('read_id', 'barcode_name', 'sequence', 'id')

        last_read = fastq_df_barcode.id.max()
    read_count = len(fastq_df_barcode)
    return fastq_df_barcode,last_read,read_count,fastasm#,fastqs_list


def custom_minimap2(newreads,previousreads,paf,firsttime=False,cores=4):
    """
    This function will run an iterative all vs all minimap2 assembly by mapping new reads to old reads, old to new and new to new.
    :param newreads: path to file containing new reads
    :param previousreads: path to file containing reads from previous iterations
    :param paf: path to file containing paf data for assembly
    :param firsttime: if True, this is the first assembly of the data
    :param cores: number of cores to use for minimap assembly
    :return: Nothing
    """
    cmd = 'minimap2 -x ava-ont -t%s %s %s >> %s ' % (cores,newreads,newreads,paf)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    status = proc.wait()
    if not firsttime:
        cmd = 'minimap2 -k15 -m100 -g10000 -r2000 --max-chain-skip 25 -t%s %s %s >> %s' % (cores, newreads, previousreads, paf)
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        status = proc.wait()
        cmd = 'minimap2 -k15 -m100 -g10000 -r2000 --max-chain-skip 25 -t%s %s %s  >> %s' % (cores, previousreads, newreads, paf)
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        status = proc.wait()

    return



@task()
def run_minimap2_assembly(job_master_id):

    job_master = JobMaster.objects.get(pk=job_master_id)

    flowcell = Flowcell.objects.get(pk=job_master.flowcell.id)

    runs = flowcell.runs.all()

    chunk_size = 25000

    fastq_df_barcode, last_read, read_count, fastasm = callfetchreads_assem(runs, chunk_size, job_master.last_read)


    if read_count > 0:
        tmp = job_master.tempfile_name

        if tmp == None:
            tmp = tempfile.NamedTemporaryFile(delete=False).name
            later = datetime.utcnow() + timedelta(days=1)
            clean_up_assembly_files.apply_async((flowcell.id,job_master_id,tmp), eta=later)

        #print (tmp)

        fastqdict = dict()

        newfastqs = 0

        for fastq in fastasm:
            if fastq.barcode not in fastqdict:
                fastqdict[fastq.barcode] = dict()
            if fastq.type not in fastqdict[fastq.barcode]:
                fastqdict[fastq.barcode][fastq.type] = []
            fastqdict[fastq.barcode][fastq.type].append([fastq.read_id, fastq.sequence])
            newfastqs += 1

        for bar in fastqdict:
            for ty in fastqdict[bar]:
                tmpfilename = tmp + bar.name + ty.name + ".gz"
                tmpfilename = tmpfilename.replace(" ", "_")
                tmpfilenameprevreads = tmp + "_prev_" + bar.name + ty.name + ".gz"
                tmpfilenameprevreads = tmpfilenameprevreads.replace(" ", "_")
                tmpfilenamepaf = tmp + bar.name + ty.name + ".paf"
                tmpfilenamepaf =tmpfilenamepaf.replace(" ", "_")
                filecontent = ""
                for fq in fastqdict[bar][ty]:
                    filecontent+=('>{}\n{}\n'.format(fq[0], fq[1]))
                #with gzip.open(tmpfilename, 'a',9) as f:
                with gzip.open(tmpfilename, 'w',9) as f:
                    f.write(filecontent.encode('utf-8'))
                    #cmd = 'minimap2 -x ava-ont -t6 %s %s | miniasm -f %s - ' % (tmpfilename, tmpfilename, tmpfilename)
                    #proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
                    #(out, err) = proc.communicate()
                    #status = proc.wait()

                if job_master.last_read==0:
                    firsttime=True
                else:
                    firsttime=False

                custom_minimap2(tmpfilename, tmpfilenameprevreads, tmpfilenamepaf, firsttime=firsttime, cores=4)

                cmd = 'cat %s >> %s' % (tmpfilename, tmpfilenameprevreads)

                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                status = proc.wait()
                os.remove(tmpfilename)

                cmd = 'miniasm -f %s %s ' % (tmpfilenameprevreads, tmpfilenamepaf)
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
                status = proc.wait()
                gfa = out.decode("utf-8")
                gfadata = gfa.splitlines()
                seqlens = []
                gfaall = ""

                for line in gfadata:
                    gfaall += line
                    line = line.strip('\n')
                    record = line.split('\t')
                    #print (line)
                    if record[0] == 'S':
                        seqlens.append(len(record[2]))

                newgfastore = GfaStore(flowcell=flowcell, barcode=bar, readtype=ty)
                newgfastore.nreads = job_master.read_count + read_count
                newgfastore.gfaformat = gfaall  # string  The whole GFA file
                newgfastore.save()
                newgfa = GfaSummary(flowcell=flowcell, barcode=bar, readtype=ty)
                newgfa.nreads = job_master.read_count + read_count
                if len(seqlens) > 0:
                    nparray = np.array(seqlens)
                    newgfa.ncontigs = len(seqlens)# int Number of contigs
                    newgfa.maxlen = max(seqlens)# int Maximum contig length
                    newgfa.minlen = min(seqlens)  # int Mininum contig length
                    newgfa.totlen = sum(seqlens)  # int Total contig length
                    newgfa.n50len = getn50(seqlens)  # int Contig N50
                    newgfa.meanlen = sum(seqlens) / len(seqlens)  # int Mean contig length
                    newgfa.allcontigs = "[%d, %d, %d, %d, %d]" % (
                    min(seqlens), np.percentile(nparray, 25), np.percentile(nparray, 50), np.percentile(nparray, 75),
                    max(seqlens))
                else:
                    newgfa.ncontigs = 0  # int Number of contigs
                    newgfa.maxlen = 0  # int Maximum contig length
                    newgfa.minlen = 0  # int Mininum contig length
                    newgfa.totlen = 0  # int Total contig length
                    newgfa.n50len = 0  # int Contig N50
                    newgfa.meanlen = 0  # int Mean contig length
                    newgfa.allcontigs = "[0,0,0,0,0]"
                newgfa.save()

    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = False
    job_master.last_read = last_read
    job_master.tempfile_name=tmp
    job_master.read_count = job_master.read_count + read_count
    job_master.save()




@task()
def run_minimap_assembly(runid, id, tmp, last_read, read_count,inputtype):

    job_master = JobMaster.objects.get(pk=id)
    # job_master.running = True
    # job_master.save()

    if last_read is None:
        last_read = 0

    if tmp == None:
        tmp = tempfile.NamedTemporaryFile(delete=False).name
        later = datetime.utcnow() + timedelta(days=1)
        clean_up_assembly_files.apply_async((runid,id,tmp), eta=later)


    logger.debug("tempfile {}".format(tmp))

    runidset = set()
    if inputtype == "flowcell":
        realflowcell = Flowcell.objects.get(pk=runid)
        flowcell_runs = realflowcell.runs.all()
        for flowcell_run in flowcell_runs:
            runidset.add(flowcell_run.id)
            # print (flowcell_run.run_id)
            # we need to get the runids that make up this run
    else:
        runidset.add(runid)

    fastqs = FastqRead.objects.filter(run_id__id__in=runidset, id__gt=int(last_read))[:10000]

    fastqdict=dict()

    newfastqs = 0

    for fastq in fastqs:
        if fastq.barcode not in fastqdict:
            fastqdict[fastq.barcode] = dict()
        if fastq.type not in fastqdict[fastq.barcode]:
            fastqdict[fastq.barcode][fastq.type] = []
        fastqdict[fastq.barcode][fastq.type].append([fastq.read_id, fastq.sequence])
        newfastqs += 1
        #read = read + '>{} \r\n{}\r\n'.format(fastq.read_id, fastq.sequence)
        #fastqdict[fastq.read_id]=fastq
        #fastqtypedict[fastq.read_id]=fastq.type
        #outtemp.write('>{}\n{}\n'.format(fastq.read_id, fastq.sequence))
        last_read = fastq.id

    #outtemp.close()
    logger.info("Added {} new fastqs".format(newfastqs))

    if newfastqs < 1000:
        JobMaster.objects.filter(pk=id).update(running=False, tempfile_name=tmp)

    else:

        totreads = read_count + newfastqs
        #print(fastqdict)
        for bar in fastqdict:
            #print(bar.name)
            for ty in fastqdict[bar]:
                #print(ty.name)
                tmpfilename = tmp+bar.name+ty.name
                tmpfilename = tmpfilename.replace(" ","_")
                outtemp = open(tmpfilename, 'a')
                for fq in fastqdict[bar][ty]:
                    outtemp.write('>{}\n{}\n'.format(fq[0],fq[1]))
                outtemp.close()
                totfq = subprocess.check_output('grep -c ">" '+tmpfilename, shell=True)
                totfq.rstrip()
                cmd = 'minimap -Sw5 -L100 -m0 -t4 %s %s | miniasm -f %s - ' % (tmpfilename, tmpfilename, tmpfilename)
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
                status = proc.wait()
                gfa = out.decode("utf-8")
                gfadata = gfa.splitlines()
                print(gfadata)
                if inputtype == "flowcell":
                    instance = Flowcell.objects.get(pk=runid)
                else:
                    instance = Run.objects.get(pk=runid)
                seqlens = []
                gfaall = ""
                for line in gfadata:
                    gfaall += line
                    line = line.strip('\n')
                    record = line.split('\t')
                    if record[0] == 'S':
                        seqlens.append(len(record[2]))
                if inputtype == "flowcell":
                    newgfastore = GfaStore(flowcell=instance, barcode=bar, readtype=ty)
                else:
                    newgfastore = GfaStore(run=instance, barcode = bar, readtype = ty)
                newgfastore.nreads = totfq
                newgfastore.gfaformat = gfaall  #   string  The whole GFA file
                newgfastore.save()
                if inputtype == "flowcell":
                    newgfa = GfaSummary(flowcell=instance, barcode = bar, readtype = ty)
                else:
                    newgfa = GfaSummary(run=instance, barcode = bar, readtype = ty)
                newgfa.nreads = totfq
                if len(seqlens) > 0:
                    nparray = np.array(seqlens)
                    newgfa.ncontigs = len(seqlens)# int Number of contigs
                    newgfa.maxlen = max(seqlens)  # int Maximum contig length
                    newgfa.minlen  = min(seqlens) # int Mininum contig length
                    newgfa.totlen = sum(seqlens)  # int Total contig length
                    newgfa.n50len = getn50(seqlens) # int Contig N50
                    newgfa.meanlen = sum(seqlens)/len(seqlens) #   int Mean contig length
                    newgfa.allcontigs = "[%d, %d, %d, %d, %d]" % (min(seqlens), np.percentile(nparray,25), np.percentile(nparray,50), np.percentile(nparray,75), max(seqlens))
                else:
                    newgfa.ncontigs = 0 # int Number of contigs
                    newgfa.maxlen = 0   # int Maximum contig length
                    newgfa.minlen  = 0  # int Mininum contig length
                    newgfa.totlen = 0   # int Total contig length
                    newgfa.n50len = 0   #   int Contig N50
                    newgfa.meanlen = 0  #   int Mean contig length
                    newgfa.allcontigs = "[0,0,0,0,0]"
                newgfa.save()


        JobMaster.objects.filter(pk=id).update(running=False, last_read=last_read, tempfile_name=tmp, read_count=totreads)


@task
def updateReadNamesOnRedis():
    print('>>> running updateReadNamesOnRedis')
    r = redis.StrictRedis(host='localhost', port=6379, db=0)

    runs = Run.objects.all()

    for run in runs:
        print('>>> run: {}'.format(run.id))

        reads = FastqRead.objects.filter(run_id=run).order_by('id')

        paginator = Paginator(reads, 25)

        key = 'run.{}.reads.number_pages'.format(run.id)

        r.set(key, paginator.num_pages)

        for page in range(1, paginator.num_pages + 1):
            print('>>> run {}, page {} of {}'.format(run.id, page, paginator.num_pages))

            key = 'run.{}.reads.page.{}'.format(run.id, page)

            result = paginator.page(page)

            result2 = set()

            for key in result:
                result2.add(key.read_id)

            r.set(key, json.dumps(list(result2)))


@task
def delete_runs():
    Run.objects.filter(to_delete=True).delete()


@task
def send_messages():
    new_messages = Message.objects.filter(delivered_date=None)

    for new_message in new_messages:

        # print('Sending message: {}'.format(new_message))

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
            # status = '@{} {}'.format(new_message.recipient.extendedopts.twitterhandle,new_message.title)
            # t.statuses.update(
            #    status=status
            # )

            message_sent = True

        if message_sent:
            print('inside message_sent')
            new_message.delivered_date = utcnow()
            new_message.save()

@task
def update_run_start_time():
    """
    This method update the field start_time of run based on the live data or
    the header of the first fastq read
    """

    runs = Run.objects.all()

    for run in runs:

        if len(run.RunDetails.all()):

            run.start_time = run.RunDetails.last().minKNOW_start_time
            origin = 'Live data'
            run.save()
            print('Updating start_time for run {} from {}'.format(run.runid, origin))
        '''
        else:
            ## This query is really slow - so - gonna move it out of here!
            #fastq = FastqRead.objects.filter(run=run).order_by('start_time').first()
            #run.start_time = fastq.start_time
            starttime = run.start_time
            if starttime is None:
                #starttime = None
                fastq = FastqRead.objects.filter(run=run)
            else:
                fastq = FastqRead.objects.filter(run=run).filter(start_time__lte=starttime)
            run.start_time = fastq.aggregate(Min('start_time'))['start_time__min']

            origin = 'Basecalled data'

        run.save()
        print('Updating start_time for run {} from {}'.format(run.runid, origin))
        '''


@task
def update_flowcell_list_details():
    """
    
    :return:
    """
    flowcell_list = Flowcell.objects.filter(is_active=True)

    for flowcell in flowcell_list:

        flowcell_job_list = JobMaster.objects.filter(flowcell=flowcell)

        for flowcell_job in flowcell_job_list:

            if flowcell_job.job_type.name == "UpdateFlowcellDetails":

                if not flowcell_job.running:

                    flowcell_job.running = True
                    flowcell_job.save()
                    update_flowcell_details(flowcell.id, flowcell_job.id)


@task
def update_flowcell_details(job_master_id):
    """
    This task updates the flowcell details (number of runs, number of reads, sample name)
    using the MinIONRunStatus records if they are available, otherwise reading from the
    fastq file summaries
    """

    """
    Refactoring this task such that:
    1. it only runs when new reads have appeared.
    2. it only processes reads since the last one was processed.
    """

    job_master = JobMaster.objects.get(pk=job_master_id)
    # job_master.running = True
    # job_master.save()

    flowcell = job_master.flowcell

    logger.info('Flowcell id: {} - Updating details of flowcell {}'.format(flowcell.id, flowcell.name))

    #
    # Get the first MinIONRunStatus for a particular flowcell
    #

    ## Seems fast enough.
    minion_run_status_first = MinIONRunStatus.objects.filter(run_id__flowcell=flowcell).order_by('minKNOW_start_time')\
        .first()

    #
    # If the MinIONRunStatus exists, than update start time and sample name
    #
    if minion_run_status_first:

        logger.info('Flowcell id: {} - There is at least one MinIONRunStatus'.format(flowcell.id))

        flowcell.start_time = minion_run_status_first.minKNOW_start_time
        flowcell.sample_name = minion_run_status_first.minKNOW_sample_name

        logger.info('Flowcell id: {} - Setting start_time to {}'.format(flowcell.id, flowcell.start_time))
        logger.info('Flowcell id: {} - Setting sample_name to {}'.format(flowcell.id, flowcell.sample_name))

    #
    # Get number of fastqreads
    #
    number_reads = 0

    for run in flowcell.runs.all():
        number_reads = number_reads + run.reads.all().count()
        if run.has_fastq:
            flowcell.has_fastq=True

    #
    # Get the job_master chancalc for this flowcell
    #
    job_master_list = JobMaster.objects.filter(flowcell=flowcell, job_type__name='Chancalc')

    number_reads_processed = flowcell.number_reads_processed

    if job_master_list.count() > 0:

        number_reads_processed = job_master_list[0].read_count

    #
    # Get the FlowcellSummaryBarcodes for a particular flowcell and for barcode_name "All reads"
    #
    flowcell_summary_list = FlowcellSummaryBarcode.objects.filter(flowcell=flowcell).filter(barcode_name='All reads')

    average_read_length = 0
    total_read_length = 0
    # number_reads_processed = 0

    logger.info('Flowcell id: {} - There is/are {} FlowcellSummaryBarcode records'.format(flowcell.id, len(flowcell_summary_list)))

    for flowcell_summary in flowcell_summary_list:
        total_read_length += flowcell_summary.total_length
        # number_reads_processed += flowcell_summary.read_count

    if number_reads > 0:
        average_read_length = total_read_length / number_reads



    logger.info('Flowcell id: {} - Total read length {}'.format(flowcell.id, total_read_length))
    logger.info('Flowcell id: {} - Number reads {}'.format(flowcell.id, number_reads))
    logger.info('Flowcell id: {} - Number reads processed {}'.format(flowcell.id, number_reads_processed))
    logger.info('Flowcell id: {} - Average read length {}'.format(flowcell.id, average_read_length))

    flowcell.average_read_length = average_read_length
    flowcell.total_read_length = total_read_length
    flowcell.number_reads = number_reads
    flowcell.number_reads_processed = number_reads_processed

    flowcell.number_runs = len(flowcell.runs.all())
    ##Faster way of doing this by using the barcode table!
    #barcode_count=0
    #runs = flowcell.runs.all()
    #for run in runs:
        #barcode_count += len(FastqRead.objects.filter(run=run).values('barcode_name').distinct())
    barcode_count = len(Barcode.objects.filter(run_id__flowcell=flowcell).values('name').distinct())-1
    flowcell.number_barcodes = barcode_count
    flowcell.save()

    logger.info('Flowcell id: {} - Number runs {}'.format(flowcell.id, flowcell.number_runs))
    logger.info('Flowcell id: {} - Number barcodes {}'.format(flowcell.id, flowcell.number_barcodes))

    #
    # Update flowcell size
    #
    ### This is a very slow query.
    ### To solve we are going to track how many reads we have looked at in the job_master - which here is the update_flow_cell task.

    runs = flowcell.runs.all()

    for run in runs:
        max_channel = FastqRead.objects.filter(run=run, id__gt=int(job_master.last_read)).aggregate(result=Max('channel'),
                                                                                 last_read=Max('id'))
        if max_channel['result'] != None and max_channel['result'] > 0:
            # we have data - so go off and use it.
            break

    # Now we need to compare the previous result with this result. To do this we need to save something else on the flowcell.

    if max_channel['result'] != None and max_channel['result'] > flowcell.max_channel:
        flowcell.max_channel = max_channel['result']
        if max_channel['result']:

            if max_channel['result'] > 512:

                flowcell.size = 3000

            elif max_channel['result'] > 126:

                flowcell.size = 512

            else:

                flowcell.size = 126

        else:

            flowcell.size = 512

    flowcell.save()

    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = False
    if max_channel['last_read'] != None:
        job_master.last_read = max_channel['last_read']
    job_master.save()

