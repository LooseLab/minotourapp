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
from django.core.mail import send_mail
from django.core.paginator import Paginator
from django.db.models import F
from django_mailgun import MailgunAPIError
from twitter import *

from alignment.models import PafRoughCov, PafStore, PafSummaryCov
from assembly.models import GfaStore, GfaSummary
from communication.utils import *
from devices.models import Flowcell
from jobs.models import JobMaster
from reads.models import Barcode, FastqRead, Run, HistogramSummary
from reads.services import (save_flowcell_histogram_summary, save_flowcell_channel_summary, save_flowcell_summary_barcode,
                            save_flowcell_statistic_barcode)
from reference.models import ReferenceInfo

from centrifuge.centrifuge import Centrifuger
import sys

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

            if flowcell_job.job_type.name == "Minimap2":

                print("trying to run alignment for flowcell {} {} {} {}".format(
                    flowcell.id,
                    flowcell_job.id,
                    flowcell_job.reference.id,
                    flowcell_job.last_read
                ))

                run_minimap2_alignment.delay(
                    flowcell.id,
                    flowcell_job.id,
                    flowcell_job.reference.id,
                    flowcell_job.last_read,
                    "flowcell"
                )

            if flowcell_job.job_type.name == "ChanCalc":

                processreads.delay(flowcell.id, flowcell_job.id, flowcell_job.last_read)

            if flowcell_job.job_type.name == "Assembly":

                inputtype = "flowcell"

                run_minimap_assembly.delay(flowcell.id, flowcell_job.id, flowcell_job.tempfile_name, flowcell_job.last_read,
                                           flowcell_job.read_count, inputtype)
            if flowcell_job.job_type.name == "Centrifuge":
                """
                    Run the Centrifuger class on the metagenomics data
                """
                print("trying to run classification for flowcell {} {} {} ".format(
                    flowcell.id,
                    flowcell_job.id,
                    flowcell_job.last_read
                ))
                print("starting centrifuge task")
                """
                starts the centrifuge instance as a celery task which should return asynchronously
                :param request:
                :return:
                """
                try:
                    run_centrifuge.delay(flowcell.id, flowcell_job.id)

                except Exception as e:
                    e = sys.exc_info()
                    print(e)

@task()
def run_centrifuge(flowcell_id, flowcell_job_id):
    print("trying centrifuge task")
    c = Centrifuger(flowcell_id, flowcell_job_id)
    c.run_centrifuge()
    print("finished task")

@task()
def processreads(flowcell_id, job_master_id, last_read):

    print('Running processreads - flowcell: {}, last read: {}, job master id: {}'.format(flowcell_id, last_read, job_master_id))

    # run = Run.objects.get(pk=runid)

    # barcode_allreads = None
    #
    # for barcode in run.barcodes.all():
    #
    #     if barcode.name == 'All reads':
    #
    #         barcode_allreads = barcode

    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = True
    job_master.save()

    fastqs = FastqRead.objects.filter(run__flowcell_id=flowcell_id).filter(id__gt=int(last_read))[:2000]

    print('found {} reads'.format(len(fastqs)))

    if len(fastqs) > 0:

        fastq_df_barcode = pd.DataFrame.from_records(fastqs.values('id', 'start_time', 'barcode__name', 'type__name', 'is_pass', 'sequence_length', 'quality_average', 'channel'))
        fastq_df_barcode['status'] = np.where(fastq_df_barcode['is_pass'] == False, 'Fail', 'Pass')
        fastq_df_barcode['start_time_truncate'] = np.array(fastq_df_barcode['start_time'], dtype='datetime64[m]')

        fastq_df_allreads = fastq_df_barcode.copy()
        fastq_df_allreads['barcode__name'] = 'All reads'

        fastq_df = fastq_df_barcode.append(fastq_df_allreads)

        #
        # Calculates statistics for RunSummaryBarcode
        #
        fastq_df_result = fastq_df.groupby(['barcode__name', 'type__name', 'is_pass']).agg({'sequence_length': ['min', 'max', 'sum', 'count'], 'quality_average': ['sum'], 'channel': ['unique']})

        fastq_df_result.reset_index().apply(lambda row: save_flowcell_summary_barcode(flowcell_id, row), axis=1)

        # fastq_df['start_time']=fastq_df['start_time'].values.astype('<M8[m]')
        #
        # Calculates statistics for RunStatisticsBarcode
        #
        fastq_df_result = fastq_df.groupby(['start_time_truncate', 'barcode__name', 'type__name', 'is_pass']).agg(
            {'sequence_length': ['min', 'max', 'sum', 'count'], 'quality_average': ['sum'], 'channel': ['unique']})

        fastq_df_result.reset_index().apply(lambda row: save_flowcell_statistic_barcode(flowcell_id, row), axis=1)

        #
        # Calculates statistics for HistogramSummary
        #
        fastq_df['bin_index'] = (fastq_df['sequence_length'] - fastq_df['sequence_length'] % HistogramSummary.BIN_WIDTH) / HistogramSummary.BIN_WIDTH

        fastq_df_result = fastq_df.groupby(['barcode__name', 'type__name', 'is_pass', 'bin_index']).agg({'sequence_length': ['sum', 'count']})

        fastq_df_result.reset_index().apply(lambda row: save_flowcell_histogram_summary(flowcell_id, row), axis=1)

        #
        # Calculates statistics for ChannelSummary
        #
        fastq_df_result = fastq_df.groupby(['channel']).agg({'sequence_length': ['sum', 'count']})

        fastq_df_result.reset_index().apply(lambda row: save_flowcell_channel_summary(flowcell_id, row), axis=1)

        last_read = fastq_df_barcode['id'].max()

    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = False
    job_master.last_read = last_read
    job_master.save()


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


def get_runidset(runid, inputtype):

    run_set = set()

    if inputtype == "flowcell":

        flowcell = Flowcell.objects.get(pk=runid)

        for run in flowcell.runs.all():

            run_set.add(run.id)

    else:

        run_set.add(runid)

    return run_set


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
        flowcell_runs = realflowcell.runs.all()
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
def run_minimap2_alignment(flowcell_id, job_master_id, reference_info_id, last_read, inputtype):

    logger.info("---> Task run_minimap2_alignment")
    logger.info("---> Parameters")
    logger.info("---> flowcell_id: {}".format(flowcell_id))
    logger.info("---> job_master_id: {}".format(job_master_id))
    logger.info("---> reference: {}".format(reference_info_id))
    logger.info("---> last_read: {}".format(last_read))
    logger.info("---> inputtype: {}".format(inputtype))

    if last_read is None:
        last_read = 0

    try:

        JobMaster.objects.filter(pk=job_master_id).update(running=True)

        REFERENCELOCATION = getattr(settings, "REFERENCELOCATION", None)

        reference = ReferenceInfo.objects.get(pk=reference_info_id)
        reference_info = ReferenceInfo.objects.get(pk=reference_info_id)

        chromdict = dict()

        chromosomes = reference.referencelines.all()

        for chromosome in chromosomes:

            chromdict[chromosome.line_name] = chromosome

        minimap2 = reference.minimap2_index_file_location

        minimap2_ref = os.path.join(REFERENCELOCATION, minimap2)

        fastqs = FastqRead.objects.filter(run__flowcell_id=flowcell_id, id__gt=int(last_read))[:1000]

        print("found fastqs: {}".format(len(fastqs)))

        read = ''
        fastqdict = dict()
        fastqtypedict = dict()
        fastqbarcodegroup = dict()
        fastqbarcode=dict()

        for fastq in fastqs:
            read = read + '>{} \r\n{}\r\n'.format(fastq.read_id, fastq.fastqreadextra.sequence)
            fastqdict[fastq.read_id] = fastq
            fastqtypedict[fastq.read_id] = fastq.type
            fastqbarcodegroup[fastq.read_id] = fastq.barcode.barcodegroup
            fastqbarcode[fastq.read_id] = fastq.barcode
            last_read = fastq.id

        minimap2_executable = getattr(settings, "MINIMAP2", None)

        if not minimap2_executable:

            print('Can not find minimap2 executable - stopping task.')


        cmd = '{} -x map-ont -t 4 --secondary=no {} -'.format(minimap2_executable, minimap2_ref)

        print(cmd)

        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                stdin=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate(input=read.encode("utf-8"))
        status = proc.wait()
        paf = out.decode("utf-8")

        pafdata = paf.splitlines()

        print(pafdata)

        #grouprun = GroupRun.objects.get(pk=runid)
        flowcell = Flowcell.objects.get(pk=flowcell_id)

        resultstore = dict()

        ### OK - we need to fix this for getting the group barcodes and not the individual barcodes.

        bulk_paf = []
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
                newpaf = PafStore(flowcell=flowcell, read=readid, read_type=typeid)
                newpafstart = PafRoughCov(flowcell=flowcell,read_type=typeid,barcode=fastqbarcode[record[0]])
                newpafend = PafRoughCov(flowcell=flowcell,read_type=typeid,barcode=fastqbarcode[record[0]])
            elif inputtype == "run":
                newpaf = PafStore(run=run, read=readid, read_type=typeid)

            newpaf.reference = reference_info

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

            newpafstart.reference = reference_info
            newpafend.reference = reference_info
            newpafstart.chromosome = chromdict[record[5]]
            newpafend.chromosome = chromdict[record[5]]
            newpafstart.p = int(record[7])
            newpafend.p = int(record[8])
            newpafstart.i = 1
            newpafend.i = -1

            bulk_paf_rough.append(newpafstart)
            bulk_paf_rough.append(newpafend)
            #newpaf.save()
            bulk_paf.append(newpaf)

            if reference_info not in resultstore:
                resultstore[reference_info] = dict()

            if chromdict[record[5]] not in resultstore[reference_info]:
                resultstore[reference_info][chromdict[record[5]]] = dict()

            #readidbarcodegroup = readid.barcode.barcodegroup
            readidbarcodegroup = fastqbarcodegroup[record[0]]

            if readidbarcodegroup not in resultstore[reference_info][chromdict[record[5]]]:
                resultstore[reference_info][chromdict[record[5]]][readidbarcodegroup] = dict()

            if typeid not in resultstore[reference_info][chromdict[record[5]]][readidbarcodegroup]:
                resultstore[reference_info][chromdict[record[5]]][readidbarcodegroup][typeid] = dict()

            if 'read' not in resultstore[reference_info][chromdict[record[5]]][readidbarcodegroup][typeid]:
                resultstore[reference_info][chromdict[record[5]]][readidbarcodegroup][typeid]['read'] = set()
                resultstore[reference_info][chromdict[record[5]]][readidbarcodegroup][typeid]['length'] = 0

            resultstore[reference_info][chromdict[record[5]]][readidbarcodegroup][typeid]['read'].add(record[0])
            resultstore[reference_info][chromdict[record[5]]][readidbarcodegroup][typeid]['length'] += int(record[3]) - int(
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
                        # if inputtype == "run":
                        #     summarycov, created2 = PafSummaryCov.objects.update_or_create(
                        #         run=runinstance,
                        #         read_type=ty,
                        #         barcode=bc,
                        #         reference=ref,
                        #         chromosome=ch,
                        #     )
                        # elif inputtype == "flowcell":
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


        jobdone = time.time()
        #print('!!!!!!!It took {} to process the resultstore.!!!!!!!!!'.format((jobdone - donepafproc)))

    except Exception as exception: #Todo: This method of detecting an error will potentially lead to duplication of analysis
        print('An error occurred when running this task.')
        print(exception)
        JobMaster.objects.filter(pk=job_master_id).update(running=False)

    #else:
    JobMaster.objects.filter(pk=job_master_id).update(running=False, last_read=last_read, read_count=F('read_count') + len(fastqs))


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

        else:

            fastq = FastqRead.objects.filter(run=run).order_by('start_time').first()
            run.start_time = fastq.start_time

            origin = 'Basecalled data'

        run.save()
        print('Updating start_time for run {} from {}'.format(run.runid, origin))
