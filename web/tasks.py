from __future__ import absolute_import, unicode_literals

import json
import os
import shutil
import subprocess
import tempfile
import time
import zipfile
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import pytz
import redis
from celery import task
from celery.utils.log import get_task_logger
from django.conf import settings
from django.core.cache import cache
from django.core.mail import send_mail
from django.core.paginator import Paginator
from django.db.models import F
from django_mailgun import MailgunAPIError
from twitter import *

from alignment.models import PafRoughCov, PafStore, PafStore_transcriptome, PafSummaryCov, PafSummaryCov_transcriptome, \
    SamStore
from assembly.models import GfaStore, GfaSummary
from communication.utils import *
from devices.models import Flowcell
from jobs.models import JobMaster, JobType
from minikraken.models import MiniKraken, ParsedKraken
from reads.models import Barcode, FastqRead, FastqReadType, FlowCellRun, Run, HistogramSummary, GroupRun
from reads.services import (save_channelsummary, save_histogramsummary,
                            save_runstatisticbarcode, save_runsummarybarcode)
from reference.models import ReferenceInfo

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


# def group_channel_presence(series):
#
#     num_rows = len(series)
#
#     if len(series) <= 0:
#         return None
#
#     num_columns = len(series[0])
#
#     new_serie = list('0' * num_columns)
#
#     for i in range(num_rows):
#
#         row = list(series[i])
#
#         for j in range(num_columns):
#
#             if row[j] == '1':
#                 new_serie[j] = '1'
#
#                 continue
#
#     return ''.join(new_serie)
#
#     # positions = [pos for pos, char in enumerate(serie) if char == '1']


@task()
def run_monitor():

    logger.info('--------------------------------')
    logger.info('Running run_monitor celery task.')
    logger.info('--------------------------------')

    grouprun_list = GroupRun.objects.all() # TODO filter only the activies

    for grouprun in grouprun_list:

        grouprun_job_list = JobMaster.objects.filter(grouprun=grouprun).filter(running=False)

        for grouprun_job in grouprun_job_list:

            # if grouprun_job.job_type.name == "Kraken":
            #
            #     run_kraken.delay(grouprun.id, grouprun_job.id, grouprun_job.last_read, "flowcell")

            if grouprun_job.job_type.name == "Minimap2":

                print("trying to run alignment for flowcell {} {} {} {}".format(
                    grouprun.id,
                    grouprun_job.id,
                    grouprun_job.reference.id,
                    grouprun_job.last_read
                ))

                run_minimap2_alignment.delay(
                    grouprun.id,
                    grouprun_job.id,
                    grouprun_job.reference.id,
                    grouprun_job.last_read,
                    "flowcell"
                )

    minion_runs = Run.objects.all() # filter(active=True)

    for minion_run in minion_runs:

        logger.debug("found a run", minion_run)

        run_jobs = JobMaster.objects.filter(run=minion_run).filter(running=False)

        for run_job in run_jobs:

            # if run_job.job_type.name == "Alignment":
            #
            #     print("trying to run bwa alignment {} {} {} {}".format(
            #         minion_run.id,
            #         run_job.id,
            #         run_job.reference.id,
            #         run_job.last_read
            #     ))
            #
            #     run_bwa_alignment.delay(minion_run.id, run_job.id, run_job.reference.id, run_job.last_read)

            # if run_job.job_type.name == "Minimap2_trans":
            #     print("trying to run transcription alignemnt {} {} {} {}".format(
            #         minion_run.id,
            #         run_job.id,
            #         run_job.reference.id,
            #         run_job.last_read
            #     ))
            #
            #     run_minimap2_transcriptome.delay(minion_run.id, run_job.id, run_job.reference.id, run_job.last_read)

            if run_job.job_type.name == "Minimap2":

                print("trying to run alignment {} {} {} {}".format(
                    minion_run.id,
                    run_job.id,
                    run_job.reference.id,
                    run_job.last_read
                ))

                run_minimap2_alignment.delay(minion_run.id, run_job.id, run_job.reference.id, run_job.last_read, "run")

            if run_job.job_type.name == "Kraken":

                run_kraken.delay(minion_run.id, run_job.id, run_job.last_read, "run")

            if run_job.job_type.name == "ProcAlign":

                proc_alignment.delay(minion_run.id, run_job.id, run_job.reference.id, run_job.last_read)

            if run_job.job_type.name == "ChanCalc":

                processreads.delay(minion_run.id, run_job.id, run_job.last_read)


@task()
def slow_monitor():

    logger.info('Running slow_monitor celery task.')

    print("Slow Monitor Called")

    # Check for uncompleted read export tasks.
    exportreadjob = JobType.objects.get(name="ExportReads")
    jobmastercollection = JobMaster.objects.filter(job_type=exportreadjob).filter(complete=False).filter(running=False)
    for jobmaster in jobmastercollection:
        #runid,id,tmp,last_read,inputtype
        export_reads.delay(jobmaster.flowcell_id,jobmaster.id,jobmaster.tempfile_name,jobmaster.last_read,"flowcell")


    testset = {}

    cachesiz = {}

    minion_runs = Run.objects.all()

    timediff = utcnow() - timedelta(days=1)

    active_runs = Run.objects.filter(active=True).distinct()

    ## We need an active flowcell measure.

    flowcell_runs = FlowCellRun.objects.filter(run__active=True).distinct()

    flowcells = set()
    # So - loop through flowcell_runs and create a job in each sub run.

    for flowcell_run in flowcell_runs:
        logger.debug("found a flowcell", flowcell_run)
        flowcells.add(flowcell_run.flowcell)


    for flowcell in flowcells:
        flowcell_jobs = JobMaster.objects.filter(flowcell=flowcell).filter(running=False)
        for flowcell_job in flowcell_jobs:
            logger.debug(flowcell_job.job_type.name)
            if flowcell_job.job_type.name == "Assembly":
                inputtype = "flowcell"
                run_minimap_assembly.delay(flowcell.id, flowcell_job.id, flowcell_job.tempfile_name, flowcell_job.last_read,
                                           flowcell_job.read_count, inputtype)

    for minion_run in active_runs:

        print("Found an active run!")

        run_jobs = JobMaster.objects.filter(run=minion_run).filter(running=False)

        for run_job in run_jobs:
            if run_job.job_type.name == "Assembly":
                print("Running Assembly")
                inputtype="run"
                run_minimap_assembly.delay(minion_run.id, run_job.id, run_job.tempfile_name, run_job.last_read, run_job.read_count,inputtype)

    for minion_run in minion_runs:
        if minion_run.last_entry() >= timediff or minion_run.last_read() >= timediff:
            cachesiz[str(minion_run.id)] = minion_run

    try:
        testset = cache.get('a-unique-key', {})

    except:

        print('a-unique-key not found')

    print('a-unique-key is {}'.format(testset))

    deleted, added = compare_two(testset, cachesiz)

    processrun(deleted, added)

    cache.set('a-unique-key', cachesiz)

    # TODO: We need someway of removing things from the dictionary which aren't still active - otherwise things will
    # persist for ever - so a compare to dictionaries.


def processrun(deleted, added):
    for run in added:
        runinstance = Run.objects.get(pk=run)
        jobinstance = JobType.objects.get(name="ChanCalc")
        runinstance.active = True
        runinstance.save()

        newjob, created = JobMaster.objects.get_or_create(run=runinstance, job_type=jobinstance, flowcell=runinstance.flowcell)

        if created is True:
            newjob.last_read = 0

        newjob.save()
        send_message([runinstance.owner], "New Active Run", "Minotour has seen a new run start on your account. This is called {}.".format(runinstance.name))

    for run in deleted:
        try:
            runinstance = Run.objects.get(pk=run)
            runinstance.active = False
            runinstance.save()

            jobinstance = JobType.objects.get(name="ChanCalc")
            JobMaster.objects.filter(job_type=jobinstance, run=runinstance).update(complete=True)
            send_message([runinstance.owner], "Run Finished",
                         "Minotour has seen a run finish on your account. This was called {}.".format(
                             runinstance.run_name))
        except Exception as exception:
            print (exception)


def compare_two(newset, cacheset):
    print("deleted keys", newset.keys() - cacheset.keys())

    deleted = newset.keys() - cacheset.keys()

    print("added keys", cacheset.keys() - newset.keys())

    added = cacheset.keys() - newset.keys()

    return deleted, added


@task()
def processreads(runid, id, last_read):

    print('>>>> Running processreads celery task.')
    print('running processreads with {} {} {}'.format(runid, id, last_read))

    run = Run.objects.get(pk=runid)

    barcode_allreads = None

    for barcode in run.barcodes.all():

        if barcode.name == 'All reads':

            barcode_allreads = barcode

    JobMaster.objects.filter(pk=id).update(running=True)

    fastqs = FastqRead.objects.filter(run_id=runid).filter(id__gt=int(last_read))[:2000]

    print('found {} reads'.format(len(fastqs)))

    if len(fastqs) > 0:

        fastq_df_barcode = pd.DataFrame.from_records(fastqs.values())

        fastq_df_allreads = fastq_df_barcode.copy()
        fastq_df_allreads['barcode_id'] = barcode_allreads.id

        fastq_df = fastq_df_barcode.append(fastq_df_allreads)

        #
        # Calculates statistics for RunSummaryBarcode
        #
        fastq_df_result = fastq_df.groupby(['barcode_id', 'type_id', 'is_pass']).agg({'sequence_length': ['min', 'max', 'sum', 'count'], 'quality_average': ['sum'], 'channel': ['unique']})

        fastq_df_result.reset_index().apply(lambda row: save_runsummarybarcode(runid, row), axis=1)

        #
        # Calculates statistics for RunStatisticsBarcode
        #
        fastq_df_result = fastq_df.groupby(['start_time', 'barcode_id', 'type_id', 'is_pass']).agg(
            {'sequence_length': ['min', 'max', 'sum', 'count'], 'quality_average': ['sum'], 'channel': ['unique']})

        fastq_df_result.reset_index().apply(lambda row: save_runstatisticbarcode(runid, row), axis=1)

        #
        # Calculates statistics for HistogramSummary
        #
        fastq_df['bin_index'] = (fastq_df['sequence_length'] - fastq_df['sequence_length'] % HistogramSummary.BIN_WIDTH) / HistogramSummary.BIN_WIDTH

        fastq_df_result = fastq_df.groupby(['barcode_id', 'type_id', 'is_pass', 'bin_index']).agg({'sequence_length': ['sum', 'count']})

        fastq_df_result.reset_index().apply(lambda row: save_histogramsummary(runid, row), axis=1)

        #
        # Calculates statistics for ChannelSummary
        #
        fastq_df_result = fastq_df.groupby(['channel']).agg({'sequence_length': ['sum', 'count']})

        fastq_df_result.reset_index().apply(lambda row: save_channelsummary(runid, row), axis=1)

        last_read = fastqs[len(fastqs) - 1].read

    JobMaster.objects.filter(pk=id).update(running=False, last_read=last_read)


@task()
def proc_alignment(runid, id, reference, last_read):

    JobMaster.objects.filter(pk=id).update(running=True)
    sams = SamStore.objects.filter(run_id=runid, id__gt=int(last_read))[:1000]
    # fp = tempfile.TemporaryFile()
    fp = open('workfile', 'w')
    for sam in sams:
        # logger.debug(sam.samline)
        fp.write(sam.samline)
        fp.write('\n')
        last_read = sam.id
    # logger.debug(fp.read())
    fp.close()
    # subprocess.run(["samtools", "faidx", "references/hep_ref.fasta"])
    subprocess.run(
        ["samtools", "view", "-bt", "/Volumes/BigElements/human_ref/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa.fai",
         "workfile", "-o", "workfile.bam"])
    # subprocess.run(["rm", "workfile"])
    subprocess.run(["samtools", "sort", "workfile.bam", "-o", "sort_workfile.bam"])
    # subprocess.run(["rm", "workfile.bam"])
    subprocess.run(["samtools", "index", "sort_workfile.bam"])
    stdoutdata = subprocess.getoutput(
        "pysamstats --type variation sort_workfile.bam  --fasta /Volumes/BigElements/human_ref/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa")
    # print("stdoutdata:\r\n " + stdoutdata)
    # subprocess.run(["rm sort_workfile.bam"])

    """
    Then we need to:
        samtools faidx references/hep_ref.fasta
        samtools view -bt references/hep_ref.fasta.fai workfile > workfile.bam
        samtools sort workfile.bam > sort_workfile.bam
        samtools index sort_workfile.bam
        pysamstats --type variation sort_workfile.bam  --fasta references/hep_ref.fasta

    After all that we just (!) parse the pysamstats lines into the existing reference table covering this information.
    """
    JobMaster.objects.filter(pk=id).update(running=False, last_read=last_read)


@task()
def test_task(string, reference):
    REFERENCELOCATION = getattr(settings, "REFERENCELOCATION", None)
    print(reference)
    Reference = ReferenceInfo.objects.get(pk=reference)
    print(type(Reference))
    print(string, REFERENCELOCATION)
    print(Reference)
    print(Reference.minimap2_index_file_location)
    lines = Reference.referencelines.all()
    for line in lines:
        print(line.id)

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
    return ("{}\n{}\n+\n{}\n".format(lineheader,fastq.fastqreadextra.sequence,fastq.fastqreadextra.quality))


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


def get_runidset(runid,inputtype):
    runidset = set()
    if inputtype == "flowcell":
        realflowcell = Flowcell.objects.get(pk=runid)
        flowcell_runs = FlowCellRun.objects.filter(flowcell=runid)
        for flowcell_run in flowcell_runs:
            runidset.add(flowcell_run.run_id)
            # print (flowcell_run.run_id)
            # we need to get the runids that make up this run
    else:
        runidset.add(runid)
    return (runidset)


@task()
def run_minimap_assembly(runid, id, tmp, last_read, read_count,inputtype):

    JobMaster.objects.filter(pk=id).update(running=True)

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
        flowcell_runs = FlowCellRun.objects.filter(flowcell=runid)
        for flowcell_run in flowcell_runs:
            runidset.add(flowcell_run.run_id)
            # print (flowcell_run.run_id)
            # we need to get the runids that make up this run
    else:
        runidset.add(runid)

    fastqs = FastqRead.objects.filter(run_id__id__in=runidset, id__gt=int(last_read))[:10000]

    fastqdict=dict()

    newfastqs = 0

    for fastq in fastqs:
        if fastq.barcode.barcodegroup not in fastqdict:
            fastqdict[fastq.barcode.barcodegroup] = dict()
        if fastq.type not in fastqdict[fastq.barcode.barcodegroup]:
            fastqdict[fastq.barcode.barcodegroup][fastq.type] = []

        fastqdict[fastq.barcode.barcodegroup][fastq.type].append([fastq.read_id, fastq.fastqreadextra.sequence])
        newfastqs += 1
        #read = read + '>{} \r\n{}\r\n'.format(fastq.read_id, fastq.fastqreadextra.sequence)
        #fastqdict[fastq.read_id]=fastq
        #fastqtypedict[fastq.read_id]=fastq.type
        #outtemp.write('>{}\n{}\n'.format(fastq.read_id, fastq.fastqreadextra.sequence))
        last_read = fastq.id

    #outtemp.close()
    print("Added {} new fastqs".format(newfastqs))

    if newfastqs < 1000:
        JobMaster.objects.filter(pk=id).update(running=False, tempfile_name=tmp)
    else:

        totreads = read_count + newfastqs

        for bar in fastqdict:
            print(bar.name)

            for ty in fastqdict[bar]:
                print(ty.name)

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

                #logger.debug("output {}".format(gfa))

                gfadata = gfa.splitlines()

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
                    newgfastore = GfaStore(flowcell=instance, barcodegroup=bar, readtype=ty)
                else:
                    newgfastore = GfaStore(run=instance, barcodegroup = bar, readtype = ty)
                newgfastore.nreads = totfq
                newgfastore.gfaformat = gfaall  #   string  The whole GFA file
                newgfastore.save()

                #### SAVE 0 IF ASSSEMBLY FAILS
                if inputtype == "flowcell":
                    newgfa = GfaSummary(flowcell=instance, barcodegroup = bar, readtype = ty)
                else:
                    newgfa = GfaSummary(run=instance, barcodegroup = bar, readtype = ty)
                newgfa.nreads = totfq
                if len(seqlens) > 0:
                    nparray = np.array(seqlens)
                    newgfa.ncontigs = len(seqlens)#	int	Number of contigs
                    newgfa.maxlen = max(seqlens)  #	int	Maximum contig length
                    newgfa.minlen  = min(seqlens) #	int	Mininum contig length
                    newgfa.totlen = sum(seqlens)  #	int	Total contig length
                    newgfa.n50len = getn50(seqlens) #    int Contig N50
                    newgfa.meanlen = sum(seqlens)/len(seqlens) #   int Mean contig length
                    newgfa.allcontigs = "[%d, %d, %d, %d, %d]" % (min(seqlens), np.percentile(nparray, 25), np.percentile(nparray, 50), np.percentile(nparray, 75), max(seqlens))
                else:
                    newgfa.ncontigs = 0 #	int	Number of contigs
                    newgfa.maxlen = 0   #	int	Maximum contig length
                    newgfa.minlen  = 0  #	int	Mininum contig length
                    newgfa.totlen = 0   #	int	Total contig length
                    newgfa.n50len = 0   #   int Contig N50
                    newgfa.meanlen = 0  #   int Mean contig length
                    newgfa.allcontigs = "[0,0,0,0,0]"
                newgfa.save()


        JobMaster.objects.filter(pk=id).update(running=False, last_read=last_read, tempfile_name=tmp, read_count=totreads)


@task
def run_minimap2_transcriptome(runid, id, reference, last_read):
    JobMaster.objects.filter(pk=id).update(running=True)
    REFERENCELOCATION = getattr(settings, "REFERENCELOCATION", None)
    Reference = ReferenceInfo.objects.get(pk=reference)
    chromdict = dict()
    chromosomes = Reference.referencelines.all()
    for chromosome in chromosomes:
        chromdict[chromosome.line_name] = chromosome
    minimap2 = Reference.minimap2_index_file_location
    minimap2_ref = os.path.join(REFERENCELOCATION, minimap2)
    # logger.debug("runid:{} last_read:{}".format(runid,last_read))
    fastqs = FastqRead.objects.filter(run_id__id=runid, id__gt=int(last_read))[:1000]
    # logger.debug("fastqs",fastqs)
    read = ''
    fastqdict = dict()
    fastqtypedict = dict()

    # logger.debug(len(fastqs))

    for fastq in fastqs:
        read = read + '>{} \r\n{}\r\n'.format(fastq.read_id, fastq.fastqreadextra.sequence)
        fastqdict[fastq.read_id] = fastq
        fastqtypedict[fastq.read_id] = fastq.type
        last_read = fastq.id
    # logger.debug(read)
    cmd = 'minimap2 -x map-ont -t 8 --secondary=no %s -' % (minimap2_ref)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            stdin=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate(input=read.encode("utf-8"))
    status = proc.wait()
    paf = out.decode("utf-8")

    pafdata = paf.splitlines()
    runinstance = Run.objects.get(pk=runid)

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
                    # logger.debug(ref,ch,bc,ty,len(resultstore[ref][ch][bc][ty]['read']),resultstore[ref][ch][bc][ty]['length'])
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
def run_minimap2_alignment(runid, job_master_id, reference, last_read, inputtype):

    logger.info("---> Task run_minimap2_alignment")
    logger.info("---> Parameters")
    logger.info("---> runid: {}".format(runid))
    logger.info("---> job_master_id: {}".format(job_master_id))
    logger.info("---> reference: {}".format(reference))
    logger.info("---> last_read: {}".format(last_read))
    logger.info("---> inputtype: {}".format(inputtype))

    if last_read is None:
        last_read = 0

    try:
    #if True:
        starttime = time.time()

        JobMaster.objects.filter(pk=job_master_id).update(running=True)

        REFERENCELOCATION = getattr(settings, "REFERENCELOCATION", None)

        Reference = ReferenceInfo.objects.get(pk=reference)

        chromdict = dict()

        chromosomes = Reference.referencelines.all()

        for chromosome in chromosomes:

            chromdict[chromosome.line_name] = chromosome

        minimap2 = Reference.minimap2_index_file_location

        minimap2_ref = os.path.join(REFERENCELOCATION, minimap2)

        fastqs = FastqRead.objects.filter(run__groupruns=runid, id__gt=int(last_read))[:1000]

        # logger.debug("fastqs",fastqs)
        read = ''
        fastqdict = dict()
        fastqtypedict = dict()
        fastqbarcodegroup = dict()
        fastqbarcode=dict()

        # logger.debug(len(fastqs))

        for fastq in fastqs:
            read = read + '>{} \r\n{}\r\n'.format(fastq.read_id, fastq.fastqreadextra.sequence)
            fastqdict[fastq.read_id] = fastq
            fastqtypedict[fastq.read_id] = fastq.type
            fastqbarcodegroup[fastq.read_id] = fastq.barcode.barcodegroup
            fastqbarcode[fastq.read_id] = fastq.barcode
            last_read = fastq.id

        # logger.debug(read)
        cmd = 'minimap2 -x map-ont -t 4 --secondary=no %s -' % (minimap2_ref)

        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(input=read.encode("utf-8"))
        status = proc.wait()
        paf = out.decode("utf-8")
        #logger.debug(paf)
        pafdata = paf.splitlines()


        doneminimaps = time.time()
        #print('!!!!!!!It took {} to run minimap.!!!!!!!!!'.format((doneminimaps-gotreadstime)))


        if inputtype =="run":
            runinstance = Run.objects.get(pk=runid)

        grouprun = GroupRun.objects.get(pk=runid)

        resultstore = dict()

        ### OK - we need to fix this for getting the group barcodes and not the individual barcodes.

        bulk_paf=[]
        bulk_paf_rough = []

        for line in pafdata:
            line = line.strip('\n')
            record = line.split('\t')
            #print(record)
            # readid = FastqRead.objects.get(read_id=record[0])
            readid = fastqdict[record[0]]
            typeid = fastqtypedict[record[0]]
            run = fastqdict[record[0]].run_id
            if inputtype == "flowcell":
                newpaf = PafStore(grouprun=grouprun, read=readid, read_type=typeid)
                newpafstart = PafRoughCov(grouprun=grouprun,read_type=typeid,barcode=fastqbarcode[record[0]],barcodegroup=fastqbarcodegroup[record[0]])
                newpafend = PafRoughCov(grouprun=grouprun,read_type=typeid,barcode=fastqbarcode[record[0]],barcodegroup=fastqbarcodegroup[record[0]])
            elif inputtype == "run":
                newpaf = PafStore(run=run, read=readid, read_type=typeid)
            newpaf.reference = Reference

            #logger.info("---> Before parsing paf record")
            #logger.info(record)

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


            newpafstart.reference=Reference
            newpafend.reference=Reference
            newpafstart.chromosome=chromdict[record[5]]
            newpafend.chromosome=chromdict[record[5]]
            newpafstart.p=int(record[7])
            newpafend.p=int(record[8])
            newpafstart.i=1
            newpafend.i=-1

            bulk_paf_rough.append(newpafstart)
            bulk_paf_rough.append(newpafend)
            #newpaf.save()
            bulk_paf.append(newpaf)

            if Reference not in resultstore:
                resultstore[Reference] = dict()
            if chromdict[record[5]] not in resultstore[Reference]:
                resultstore[Reference][chromdict[record[5]]] = dict()
            #readidbarcodegroup = readid.barcode.barcodegroup
            readidbarcodegroup = fastqbarcodegroup[record[0]]
            if readidbarcodegroup not in resultstore[Reference][chromdict[record[5]]]:
                resultstore[Reference][chromdict[record[5]]][readidbarcodegroup] = dict()
            if typeid not in resultstore[Reference][chromdict[record[5]]][readidbarcodegroup]:
                resultstore[Reference][chromdict[record[5]]][readidbarcodegroup][typeid] = dict()
            if 'read' not in resultstore[Reference][chromdict[record[5]]][readidbarcodegroup][typeid]:
                resultstore[Reference][chromdict[record[5]]][readidbarcodegroup][typeid]['read'] = set()
                resultstore[Reference][chromdict[record[5]]][readidbarcodegroup][typeid]['length'] = 0
            resultstore[Reference][chromdict[record[5]]][readidbarcodegroup][typeid]['read'].add(record[0])
            resultstore[Reference][chromdict[record[5]]][readidbarcodegroup][typeid]['length'] += int(record[3]) - int(
                record[2]) + 1
        PafStore.objects.bulk_create(bulk_paf)
        PafRoughCov.objects.bulk_create(bulk_paf_rough)
        donepafproc = time.time()

        #print('!!!!!!!It took {} to parse the paf.!!!!!!!!!'.format((donepafproc-doneminimaps)))

        for ref in resultstore:
            for ch in resultstore[ref]:
                for bc in resultstore[ref][ch]:
                    for ty in resultstore[ref][ch][bc]:
                        # logger.debug(ref,ch,bc,ty,len(resultstore[ref][ch][bc][ty]['read']),resultstore[ref][ch][bc][ty]['length'])
                        if inputtype == "run":
                            summarycov, created2 = PafSummaryCov.objects.update_or_create(
                                run=runinstance,
                                read_type=ty,
                                barcode=bc,
                                reference=ref,
                                chromosome=ch,
                            )
                        elif inputtype == "flowcell":
                            summarycov, created2 = PafSummaryCov.objects.update_or_create(
                                flowcell=grouprun,
                                read_type=ty,
                                barcodegroup=bc,
                                reference=ref,
                                chromosome=ch,
                            )
                        summarycov.read_count += len(resultstore[ref][ch][bc][ty]['read'])
                        summarycov.cumu_length += resultstore[ref][ch][bc][ty]['length']
                        summarycov.save()


        jobdone = time.time()
        #print('!!!!!!!It took {} to process the resultstore.!!!!!!!!!'.format((jobdone - donepafproc)))

    except Exception as exception:
        print('An error occurred when running this task.')
        print(exception)
        JobMaster.objects.filter(pk=job_master_id).update(running=False)

    #else:
    JobMaster.objects.filter(pk=job_master_id).update(running=False, last_read=last_read, read_count=F('read_count') + len(fastqs))


@task()
def run_kraken(runid, id, last_read, inputtype):
    JobMaster.objects.filter(pk=id).update(running=True)

    runidset=set()
    if inputtype == "flowcell":
        flowcell_runs = FlowCellRun.objects.filter(flowcell=runid)
        for flowcell_run in flowcell_runs:
            runidset.add(flowcell_run.run_id)
        #we need to get the runids that make up this run
    else:
        runidset.add(runid)

    fastqs = FastqRead.objects.filter(run_id__in=runidset).filter(id__gt=int(last_read))[:5000]
    krakrun = Kraken()

    print (runid, id, last_read, inputtype)

    # IMPORTANT - NEEDS FIXING FOR BARCODES AND READ TYPES
    if len(fastqs) > 0:
        read = ''
        readdict = dict()
        typedict = dict()
        barcodedict = dict()
        for fastq in fastqs:
            read = '{}>{} \r\n{}\r\n'.format(read, fastq, fastq.fastqreadextra.sequence)
            # logger.debug(fastq.id)
            readdict[fastq.read_id] = fastq.id
            last_read = fastq.id
        krakrun.write_seqs(read.encode("utf-8"))
        krakenoutput = krakrun.run().decode('utf-8').split('\n')
        if inputtype == "run":
            runinstance = Run.objects.get(pk=runid)
        elif inputtype == "flowcell":
            runinstance = Flowcell.objects.get(pk=runid)
        for line in krakenoutput:
            if len(line) > 0:
                # logger.debug(line)
                krakenbits = line.split("\t")

                # logger.debug(readdict[krakenbits[1]],krakenbits[1])
                if inputtype == "run":
                    newkrak = MiniKraken(run_id=runid, read_id=int(readdict[krakenbits[1]]))
                elif inputtype == "flowcell":
                    newkrak = MiniKraken(flowcell_id=runid, read_id=int(readdict[krakenbits[1]]))
                newkrak.krakenstatus = str(krakenbits[0])
                newkrak.krakentaxid = int(krakenbits[2])
                newkrak.krakenseqlen = int(krakenbits[3])
                newkrak.krakenlca = str(krakenbits[4])
                newkrak.save()

        # To get the data about the run - specifically which barcodes it has:
        barcodes = Barcode.objects.filter(run_id__in=runidset)

        barcodeset=set()
        for bar in barcodes:
            barcodeset.add(bar.barcodegroup)

        # To get the data about the run - specifcally read types:
        types = FastqReadType.objects.all()

        ####Up to this point should work with either flowcell or run concept.
        ####after this we have to resolve how to collapse barcode types - which is complex
        ####as the system creates multiple barcodes for the same barcode.

        for bar in barcodeset:
            for type in types:
                # logger.debug(Bar,Type)
                # MiniKrak = MiniKraken.objects.filter(run_id=runid).filter(barcode=Bar).filter(read_type=Type)
                if inputtype == "run":
                    MiniKrak = MiniKraken.objects.filter(run=runinstance).filter(read__barcode__barcodegroup=bar).filter(read__type=type)
                elif inputtype =="flowcell":
                    MiniKrak = MiniKraken.objects.filter(flowcell=runinstance).filter(
                        read__barcode__barcodegroup=bar).filter(read__type=type)
                if len(MiniKrak) > 0:
                    kraken = ""
                    for krak in MiniKrak:
                        kraken = (
                        "{}{}\t{}\t{}\t{}\t{}\n".format(kraken, krak.krakenstatus, krak.read_id, krak.krakentaxid,
                                                        krak.krakenseqlen, krak.krakenlca))
                    krakrun.write_kraken(kraken.encode("utf-8"))
                    output = krakrun.process().decode('utf-8').split('\n')
                    counter = 0
                    rankdict = dict()
                    rankdict[-1] = 'Input'
                    for line in output:
                        if len(line) > 0:
                            linebits = line.split("\t")
                            # logger.debug(linebits)
                            # logger.debug(float(linebits[0].lstrip(' ')))

                            if inputtype == "run":
                                parsekrak, created = ParsedKraken.objects.get_or_create(run_id=runid, NCBItaxid=linebits[4],
                                                                                    type=type, barcode=bar)
                            elif inputtype =="flowcell":
                                parsekrak, created = ParsedKraken.objects.get_or_create(flowcell_id=runid, NCBItaxid=linebits[4],
                                                                                    type=type, barcode=bar)
                            parsekrak.percentage = float(linebits[0].lstrip(' '))
                            parsekrak.rootreads = int(linebits[1])
                            parsekrak.directreads = int(linebits[2])
                            parsekrak.rank = str(linebits[3])
                            parsekrak.sci_name = str(linebits[5].lstrip(' '))
                            parsekrak.indentation = int((len(linebits[5]) - len(linebits[5].lstrip(' '))) / 2)
                            parsekrak.orderin = counter
                            rankdict[int((len(linebits[5]) - len(linebits[5].lstrip(' '))) / 2)] = str(
                                linebits[5].lstrip(' '))
                            parsekrak.parent = rankdict[int((len(linebits[5]) - len(linebits[5].lstrip(' '))) / 2) - 1]
                            parsekrak.save()
                            counter += 1
                            # for krak in MiniKrak:
                            # logger.debug(krak)

    krakrun.finish()
    JobMaster.objects.filter(pk=id).update(running=False, last_read=last_read, read_count=F('read_count') + len(fastqs))


class Kraken():
    """
    This class assumes that kraken is available in the command line. That isn't likely to be the case.
    Should we specify a folder that contain utilities needed by minoTour?
    We also need a folder that contains the reference databases.
    In a sense we have this with the reference collections. Lets use that!
    """
    def __init__(self):
        self.tmpfile = tempfile.NamedTemporaryFile(suffix=".fa")
        self.krakenfile = tempfile.NamedTemporaryFile(suffix=".out")
        self.REFERENCELOCATION = getattr(settings, "REFERENCELOCATION", None)
        self.krakenlocation = os.path.join(self.REFERENCELOCATION, 'minikraken_20141208')

    def write_seqs(self, seqs):
        # self.tmpfile.write("\n".join(seqs))
        self.tmpfile.write(seqs)

    def write_kraken(self, kraken):
        # print "printing the kraken"
        # print "".join(kraken)
        self.krakenfile.write(kraken)

    def read_kraken(self):
        for line in self.krakenfile:
            print(line)

    def process(self):
        print('kraken-report --db '+ self.krakenlocation + ' ' + self.krakenfile.name)
        p1 = subprocess.Popen('kraken-report --db '+ self.krakenlocation + ' ' + self.krakenfile.name,
                              shell=True, stdout=subprocess.PIPE)
        (out, err) = p1.communicate()
        return out

    def process2(self):
        print('kraken-mpa-report --db '+ self.krakenlocation + ' ' + self.krakenfile.name)
        p1 = subprocess.Popen('kraken-mpa-report --db '+ self.krakenlocation + ' ' + self.krakenfile.name,shell=True, stdout=subprocess.PIPE)
        (out, err) = p1.communicate()
        return out

    def run(self):
        print('kraken --quick --db '+ self.krakenlocation + '  --fasta-input --preload  ' + self.tmpfile.name)  # +'  #| kraken-translate --db /Volumes/SSD/kraken/minikraken_20141208/ $_'
        p1 = subprocess.Popen('kraken --quick --db '+ self.krakenlocation + ' --fasta-input --preload  ' + self.tmpfile.name + '  ',shell=True, stdout=subprocess.PIPE)
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

    # logger.debug(runid,id,reference,last_read)
    fastq_list = FastqRead.objects.filter(
        run_id=runid,
        id__gt=last_read
    )[:1000]

    referenceinfo = ReferenceInfo.objects.get(pk=referenceinfo_id)

    for fastq in fastq_list:
        # logger.debug(fastq.fastqreadextra.sequence)
        bwaindex = referenceinfo.filename
        read = '>{} \r\n{}\r\n'.format(fastq, fastq.fastqreadextra.sequence)
        cmd = 'bwa mem -x ont2d %s -' % (bwaindex)
        # logger.debug(cmd)
        # logger.debug(read)
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
                    runinstance = Run.objects.get(pk=runid)
                    newsam = SamStore(run_id=runinstance, read_id=fastq)
                    newsam.samline = line
                    newsam.save()
                    last_read = fastq.id

    job.running = False
    job.last_read = last_read
    job.save()


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
