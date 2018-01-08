from __future__ import absolute_import, unicode_literals

import json
import os
import subprocess
import tempfile
import numpy as np
from datetime import datetime, timedelta
import time

import pytz
import redis
from celery import task
from celery.utils.log import get_task_logger
from celery.signals import celeryd_init
from django.conf import settings
from django.core.cache import cache
from django.core.mail import send_mail
from django.core.paginator import Paginator
from django.db.models import F
from django_mailgun import MailgunAPIError
from twitter import *



from alignment.models import PafStore, PafStore_transcriptome, PafRoughCov
from alignment.models import PafSummaryCov, PafSummaryCov_transcriptome
from alignment.models import SamStore
from communication.models import Message
from minikraken.models import MiniKraken, ParsedKraken
from reads.models import Barcode, FlowCellRun, FlowCell
from reads.models import ChannelSummary
from reads.models import FastqRead, FastqReadExtra
from reads.models import HistogramSummary
from reads.models import JobMaster
from reads.models import JobType
from reads.models import MinIONRun, FastqReadType
from reads.models import RunStatisticBarcode
from reads.models import RunSummaryBarcode
from reference.models import ReferenceInfo
from assembly.models import GfaStore
from assembly.models import GfaSummary

from communication.utils import *

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

    flowcell_runs = FlowCellRun.objects.filter(run__active=True).distinct()
    logger.debug("!!!!!!!!!!!!!!! FLOWCELL RUNS !!!!!!!!!!!!!!!!!!!!")
    logger.debug(flowcell_runs)
    flowcells=set()
    #So - loop through flowcell_runs and create a job in each sub run.

    for flowcell_run in flowcell_runs:
        logger.debug("found a flowcell", flowcell_run)
        flowcells.add(flowcell_run.flowcell)

    for flowcell in flowcells:
        flowcell_jobs = JobMaster.objects.filter(flowcell=flowcell).filter(running=False)
        for flowcell_job in flowcell_jobs:
            logger.debug(flowcell_job.job_type.name)
            if flowcell_job.job_type.name == "Kraken":
                run_kraken.delay(flowcell.id,flowcell_job.id,flowcell_job.last_read,"flowcell")

            if flowcell_job.job_type.name == "Minimap2":
                print("trying to run alignment for flowcell {} {} {} {}".format(
                    flowcell.id,
                    flowcell_job.id,
                    flowcell_job.reference.id,
                    flowcell_job.last_read
                ))

                run_minimap2_alignment.delay(flowcell.id, flowcell_job.id, flowcell_job.reference.id, flowcell_job.last_read, "flowcell")



    minion_runs = MinIONRun.objects.filter(active=True).distinct()

    for minion_run in minion_runs:
        logger.debug("found a run", minion_run)
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
                print("trying to run transcription alignemnt {} {} {} {}".format(
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

    testset = {}

    cachesiz = {}

    minion_runs = MinIONRun.objects.all()

    timediff = utcnow() - timedelta(days=1)

    active_runs = MinIONRun.objects.filter(active=True).distinct()

    for minion_run in active_runs:

        print("Found an active run!")

        run_jobs = JobMaster.objects.filter(run=minion_run).filter(running=False)

        for run_job in run_jobs:
            if run_job.job_type.name == "Assembly":
                print("Running Assembly")
                run_minimap_assembly.delay(minion_run.id, run_job.id, run_job.tempfile_name, run_job.last_read, run_job.read_count)

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
        runinstance = MinIONRun.objects.get(pk=run)
        jobinstance = JobType.objects.get(name="ChanCalc")
        runinstance.active = True
        runinstance.save()

        newjob, created = JobMaster.objects.get_or_create(run=runinstance, job_type=jobinstance)

        if created is True:
            newjob.last_read = 0

        newjob.save()
        send_message([runinstance.owner],"New Active Run","Minotour has seen a new run start on your account. This is called {}.".format(runinstance.run_name))

    for run in deleted:
        runinstance = MinIONRun.objects.get(pk=run)
        runinstance.active = False
        runinstance.save()

        jobinstance = JobType.objects.get(name="ChanCalc")
        JobMaster.objects.filter(job_type=jobinstance, run=runinstance).update(complete=True)
        send_message([runinstance.owner], "Run Finished",
                     "Minotour has seen a run finish on your account. This was called {}.".format(
                         runinstance.run_name))


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
    JobMaster.objects.filter(pk=id).update(running=True)

    fastqs = FastqRead.objects.filter(run_id=runid, id__gt=int(last_read))[:2000]

    chanstore = dict()

    histstore = dict()

    histstore['All reads'] = {}

    sumstore = dict()
    sumstoresum = dict()

    barstore = set()

    for fastq in fastqs:
        # logger.debug(fastq)
        barstore.add(fastq.barcode)
        # print('>>>> barcode: {}'.format(fastq.barcode))
        if fastq.channel not in chanstore.keys():
            chanstore[fastq.channel] = dict()
            chanstore[fastq.channel]['count'] = 0
            chanstore[fastq.channel]['length'] = 0
        chanstore[fastq.channel]['count'] += 1
        chanstore[fastq.channel]['length'] += fastq.sequence_length

        if fastq.barcode.name not in histstore.keys():
            histstore[fastq.barcode.name] = {}

        if fastq.type.name not in histstore[fastq.barcode.name].keys():
            histstore[fastq.barcode.name][fastq.type.name] = {}

        if fastq.type.name not in histstore['All reads'].keys():
            histstore['All reads'][fastq.type.name] = {}

        bin_width = findbin(fastq.sequence_length)

        if bin_width not in histstore[fastq.barcode.name][fastq.type.name].keys():
            histstore[fastq.barcode.name][fastq.type.name][bin_width] = {}

            histstore[fastq.barcode.name][fastq.type.name][bin_width]['count'] = 0
            histstore[fastq.barcode.name][fastq.type.name][bin_width]['length'] = 0

        if bin_width not in histstore['All reads'][fastq.type.name].keys():
            histstore['All reads'][fastq.type.name][bin_width] = {}

            histstore['All reads'][fastq.type.name][bin_width]['count'] = 0
            histstore['All reads'][fastq.type.name][bin_width]['length'] = 0

        histstore[fastq.barcode.name][fastq.type.name][bin_width]['count'] += 1
        histstore[fastq.barcode.name][fastq.type.name][bin_width]['length'] += fastq.sequence_length

        histstore['All reads'][fastq.type.name][bin_width]['count'] += 1
        histstore['All reads'][fastq.type.name][bin_width]['length'] += fastq.sequence_length

        ### Adding in the current post save behaviour:
        ### This is really inefficient - we are hitting the database too often.
        ### We need to store the relevant data in a single dictionary and then process it in to the database.
        ### Key elements are the time and the barcode?

        ipn_obj = fastq

        barcode = ipn_obj.barcode

        barcode_all_reads = ipn_obj.run_id.barcodes.filter(name='All reads').first()

        tm = ipn_obj.start_time

        tm = tm - timedelta(
            minutes=(tm.minute % 1) - 1,
            seconds=tm.second,
            microseconds=tm.microsecond
        )

        ### At this point we know the read time (tm) we know the barcode and we have the object.

        if ipn_obj.run_id not in sumstore.keys():
            sumstore[ipn_obj.run_id]=dict()
            sumstoresum[ipn_obj.run_id]=dict()
            
        if ipn_obj.type not in sumstore[ipn_obj.run_id].keys():
            sumstore[ipn_obj.run_id][ipn_obj.type] = dict()
            sumstoresum[ipn_obj.run_id][ipn_obj.type] = dict()


        if barcode_all_reads not in sumstoresum[ipn_obj.run_id][ipn_obj.type].keys():
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode_all_reads] = dict()
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode_all_reads]["channels"] = '0' * 3000
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode_all_reads]["pass_length"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode_all_reads]["pass_max_length"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode_all_reads]["pass_min_length"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode_all_reads]["pass_count"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode_all_reads]["total_length"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode_all_reads]["max_length"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode_all_reads]["min_length"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode_all_reads]["read_count"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode_all_reads]["quality_sum"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode_all_reads]["pass_quality_sum"] = 0

        if barcode not in sumstoresum[ipn_obj.run_id][ipn_obj.type].keys():
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode] = dict()
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode]["channels"] = '0' * 3000
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode]["pass_length"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode]["pass_max_length"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode]["pass_min_length"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode]["pass_count"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode]["total_length"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode]["max_length"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode]["min_length"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode]["read_count"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode]["quality_sum"] = 0
            sumstoresum[ipn_obj.run_id][ipn_obj.type][barcode]["pass_quality_sum"] = 0


        if tm not in sumstore[ipn_obj.run_id][ipn_obj.type].keys():
            sumstore[ipn_obj.run_id][ipn_obj.type][tm]=dict()

        if barcode_all_reads not in sumstore[ipn_obj.run_id][ipn_obj.type][tm].keys():
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_all_reads] = dict()
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_all_reads]["channels"] = '0' * 3000
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_all_reads]["pass_length"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_all_reads]["pass_max_length"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_all_reads]["pass_min_length"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_all_reads]["pass_count"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_all_reads]["total_length"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_all_reads]["max_length"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_all_reads]["min_length"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_all_reads]["read_count"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_all_reads]["quality_sum"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_all_reads]["pass_quality_sum"] = 0

        if barcode not in sumstore[ipn_obj.run_id][ipn_obj.type][tm].keys():
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode] = dict()
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode]["channels"] = '0' * 3000
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode]["pass_length"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode]["pass_max_length"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode]["pass_min_length"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode]["pass_count"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode]["total_length"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode]["max_length"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode]["min_length"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode]["read_count"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode]["quality_sum"] = 0
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode]["pass_quality_sum"] = 0



        sumstore = update_local_dict(sumstore,ipn_obj,tm,barcode_all_reads)

        '''
        obj1, created1 = RunSummaryBarcode.objects.update_or_create(
            run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=barcode_all_reads
        )
        update_sum_stats(obj1, ipn_obj)
        '''

        sumstoresum = update_local_dict_sum(sumstoresum, ipn_obj, barcode_all_reads)
        # all reads and barcodes are saved on RunStatisticBarcode
        '''
        obj3, created3 = RunStatisticBarcode.objects.update_or_create(
            run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=barcode_all_reads, sample_time=tm
        )
        update_sum_stats(obj3, ipn_obj)
        '''

        if barcode is not None and barcode.name != '':
            sumstore = update_local_dict(sumstore, ipn_obj, tm, barcode)
            '''
            obj2, created2 = RunSummaryBarcode.objects.update_or_create(
                run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=barcode
            )
            update_sum_stats(obj2, ipn_obj)
            '''
            sumstoresum = update_local_dict_sum(sumstoresum, ipn_obj, barcode)
            '''
            obj3, created3 = RunStatisticBarcode.objects.update_or_create(
                run_id=ipn_obj.run_id, type=ipn_obj.type, barcode=barcode, sample_time=tm
            )
            update_sum_stats(obj3, ipn_obj)
            '''
        last_read = fastq.id

    runinstance = MinIONRun.objects.get(pk=runid)

    logger.debug("!!!!!!!!! SUMSTORE ¢¢¢¢", len(sumstore))
    logger.debug("!!!!!!!!! SUMSTORESUM ¢¢¢¢", len(sumstoresum))

    for run in sumstore:
        for type_ in sumstore[run]:
            for tm in sumstore[run][type_]:
                for barcode in sumstore[run][type_][tm]:
                    obj, created1 = RunStatisticBarcode.objects.update_or_create(
                        run_id=run, type=type_, barcode=barcode, sample_time=tm
                    )
                    #update_sum_stats(obj1, sumstore[run][type_][tm][barcode])
                    obj.pass_length += sumstore[run][type_][tm][barcode]["pass_length"]
                    obj.pass_quality_sum += sumstore[run][type_][tm][barcode]["pass_quality_sum"]

                    if sumstore[run][type_][tm][barcode]["pass_max_length"] > obj.pass_max_length:
                        obj.pass_max_length = sumstore[run][type_][tm][barcode]["pass_max_length"]

                    if obj.pass_min_length == 0:
                        obj.pass_min_length = sumstore[run][type_][tm][barcode]["pass_min_length"]

                    if sumstore[run][type_][tm][barcode]["pass_min_length"] < obj.pass_min_length:
                        obj.pass_min_length = sumstore[run][type_][tm][barcode]["pass_min_length"]

                    obj.pass_count += sumstore[run][type_][tm][barcode]["pass_count"]

                    obj.total_length += sumstore[run][type_][tm][barcode]["total_length"]
                    obj.quality_sum += sumstore[run][type_][tm][barcode]["quality_sum"]

                    if sumstore[run][type_][tm][barcode]["max_length"] > obj.max_length:
                        obj.max_length = sumstore[run][type_][tm][barcode]["max_length"]
    
                    if obj.min_length == 0:
                        obj.min_length = sumstore[run][type_][tm][barcode]["min_length"]
    
                    if sumstore[run][type_][tm][barcode]["min_length"] < obj.min_length:
                        obj.min_length = sumstore[run][type_][tm][barcode]["min_length"]
    
                    #
                    # channel_presence is a 512 characters length string containing 0
                    # for each channel (id between 1 - 512) seen, we set the position on
                    # the string to 1. afterwards, we can calculate the number of active
                    # channels in a particular minute.
                    #

                    channel = sumstore[run][type_][tm][barcode]["channels"]
                    channel_sequence = obj.channel_presence
                    channel_sequence_list = list(channel_sequence)
                    for i, val in enumerate(channel):
                        if int(val) == 1:
                            channel_sequence_list[i] = '1'
                    obj.channel_presence = ''.join(channel_sequence_list)

                    obj.read_count += sumstore[run][type_][tm][barcode]["read_count"]
                    obj.save()

    for run in sumstoresum:
        for type_ in sumstoresum[run]:
            for barcode in sumstoresum[run][type_]:
                obj, created1 = RunSummaryBarcode.objects.update_or_create(
                    run_id=run, type=type_, barcode=barcode
                )
                # update_sum_stats(obj1, sumstoresum[run][type_][barcode])
                obj.pass_length += sumstoresum[run][type_][barcode]["pass_length"]
                obj.pass_quality_sum += sumstoresum[run][type_][barcode]["pass_quality_sum"]

                if sumstoresum[run][type_][barcode]["pass_max_length"] > obj.pass_max_length:
                    obj.pass_max_length = sumstoresum[run][type_][barcode]["pass_max_length"]

                if obj.pass_min_length == 0:
                    obj.pass_min_length = sumstoresum[run][type_][barcode]["pass_min_length"]

                if sumstoresum[run][type_][barcode]["pass_min_length"] < obj.pass_min_length:
                    obj.pass_min_length = sumstoresum[run][type_][barcode]["pass_min_length"]

                obj.pass_count += sumstoresum[run][type_][barcode]["pass_count"]

                obj.total_length += sumstoresum[run][type_][barcode]["total_length"]
                obj.quality_sum += sumstoresum[run][type_][barcode]["quality_sum"]

                if sumstoresum[run][type_][barcode]["max_length"] > obj.max_length:
                    obj.max_length = sumstoresum[run][type_][barcode]["max_length"]

                if obj.min_length == 0:
                    obj.min_length = sumstoresum[run][type_][barcode]["min_length"]

                if sumstoresum[run][type_][barcode]["min_length"] < obj.min_length:
                    obj.min_length = sumstoresum[run][type_][barcode]["min_length"]

                #
                # channel_presence is a 512 characters length string containing 0
                # for each channel (id between 1 - 512) seen, we set the position on
                # the string to 1. afterwards, we can calculate the number of active
                # channels in a particular minute.
                #

                channel = sumstoresum[run][type_][barcode]["channels"]
                channel_sequence = obj.channel_presence
                channel_sequence_list = list(channel_sequence)
                for i,val in enumerate(channel):
                    if int(val)==1:
                        channel_sequence_list[i]='1'
                obj.channel_presence = ''.join(channel_sequence_list)

                obj.read_count += sumstoresum[run][type_][barcode]["read_count"]
                obj.save()

    for chan in chanstore:
        channel, created = ChannelSummary.objects.get_or_create(run_id=runinstance, channel_number=int(chan))
        channel.read_count += chanstore[chan]['count']
        channel.read_length += chanstore[chan]['length']
        channel.save()

    for barcodename in histstore.keys():

        barcode = Barcode.objects.filter(run=runinstance, name=barcodename).first()

        for read_type_name in histstore[barcodename].keys():
            read_type = FastqReadType.objects.get(name=read_type_name)
            # print('>>>> {}'.format(read_type_name))

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
    logger.debug("Finished this thang!")


def update_local_dict(sumstore,ipn_obj,tm,barcode_to_do):
    channel = ipn_obj.channel
    channel_sequence = sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["channels"]
    channel_sequence_list = list(channel_sequence)
    ##Temporary fix to enable promethION data upload
    #if channel <= 512:
    channel_sequence_list[channel - 1] = '1'
    sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["channels"] = ''.join(channel_sequence_list)

    if ipn_obj.is_pass:


        sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["pass_length"] += ipn_obj.sequence_length
        sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["pass_quality_sum"] += ipn_obj.quality_average

        if ipn_obj.sequence_length > sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["pass_max_length"]:
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["pass_max_length"] = ipn_obj.sequence_length

        if sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["pass_min_length"] == 0:
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["pass_min_length"] = ipn_obj.sequence_length

        if ipn_obj.sequence_length < sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["pass_min_length"]:
            sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["pass_min_length"] = ipn_obj.sequence_length

        sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["pass_count"] += 1

    sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["total_length"] += ipn_obj.sequence_length

    sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["quality_sum"] += ipn_obj.quality_average

    if ipn_obj.sequence_length > sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["max_length"]:
        sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["max_length"] = ipn_obj.sequence_length

    if sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["min_length"] == 0:
        sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["min_length"] = ipn_obj.sequence_length

    if ipn_obj.sequence_length < sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["min_length"]:
        sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["min_length"] = ipn_obj.sequence_length

    sumstore[ipn_obj.run_id][ipn_obj.type][tm][barcode_to_do]["read_count"] += 1

    return sumstore


def update_local_dict_sum(sumstore, ipn_obj, barcode_to_do):
    channel = ipn_obj.channel
    channel_sequence = sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["channels"]
    channel_sequence_list = list(channel_sequence)
    ##Temporary fix to enable promethION data upload
    #if channel <= 512:
    channel_sequence_list[channel - 1] = '1'
    sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["channels"] = ''.join(channel_sequence_list)

    if ipn_obj.is_pass:

        sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["pass_length"] += ipn_obj.sequence_length
        sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["pass_quality_sum"] += ipn_obj.quality_average

        if ipn_obj.sequence_length > sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["pass_max_length"]:
            sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["pass_max_length"] = ipn_obj.sequence_length

        if sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["pass_min_length"] == 0:
            sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["pass_min_length"] = ipn_obj.sequence_length

        if ipn_obj.sequence_length < sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["pass_min_length"]:
            sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["pass_min_length"] = ipn_obj.sequence_length

        sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["pass_count"] += 1

    sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["total_length"] += ipn_obj.sequence_length
    sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["quality_sum"] += ipn_obj.quality_average
    if ipn_obj.sequence_length > sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["max_length"]:
        sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["max_length"] = ipn_obj.sequence_length

    if sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["min_length"] == 0:
        sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["min_length"] = ipn_obj.sequence_length

    if ipn_obj.sequence_length < sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["min_length"]:
        sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["min_length"] = ipn_obj.sequence_length

    sumstore[ipn_obj.run_id][ipn_obj.type][barcode_to_do]["read_count"] += 1

    return sumstore

def update_sum_stats(obj, ipn_obj):
    if ipn_obj.is_pass:

        obj.pass_length += ipn_obj.sequence_length

        if ipn_obj.sequence_length > obj.pass_max_length:
            obj.pass_max_length = ipn_obj.sequence_length

        if obj.pass_min_length == 0:
            obj.pass_min_length = ipn_obj.sequence_length

        if ipn_obj.sequence_length < obj.pass_min_length:
            obj.pass_min_length = ipn_obj.sequence_length

        obj.pass_count += 1

    obj.total_length += ipn_obj.sequence_length

    if ipn_obj.sequence_length > obj.max_length:
        obj.max_length = ipn_obj.sequence_length

    if obj.min_length == 0:
        obj.min_length = ipn_obj.sequence_length

    if ipn_obj.sequence_length < obj.min_length:
        obj.min_length = ipn_obj.sequence_length

    #
    # channel_presence is a 512 characters length string containing 0
    # for each channel (id between 1 - 512) seen, we set the position on
    # the string to 1. afterwards, we can calculate the number of active
    # channels in a particular minute.
    #
    channel = ipn_obj.channel
    channel_sequence = obj.channel_presence
    channel_sequence_list = list(channel_sequence)
    channel_sequence_list[channel - 1] = '1'
    obj.channel_presence = ''.join(channel_sequence_list)

    obj.read_count += 1
    obj.save()


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
def run_minimap_assembly(runid, id, tmp, last_read, read_count):

    JobMaster.objects.filter(pk=id).update(running=True)
    logger.debug("hello teri {}".format(runid))
    logger.debug("tempfile {}".format(tmp))
    if tmp == None:
        tmp = tempfile.NamedTemporaryFile(delete=False).name
    logger.debug("tempfile {}".format(tmp))

    fastqs = FastqRead.objects.filter(run_id__id=runid, id__gt=int(last_read))[:10000]
    #logger.debug("fastqs",fastqs)
    #read = ''
    #fastqdict=dict()
    #fastqtypedict=dict()
    #outtemp = open(tmp, 'a')

    fastqdict=dict()

    newfastqs = 0

    for fastq in fastqs:
        if fastq.barcode not in fastqdict:
            fastqdict[fastq.barcode] = dict()
        if fastq.type not in fastqdict[fastq.barcode]:
            fastqdict[fastq.barcode][fastq.type] = []

        fastqdict[fastq.barcode][fastq.type].append([fastq.read_id, fastq.fastqreadextra.sequence])
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

                runinstance = MinIONRun.objects.get(pk=runid)

                seqlens = []
                gfaall = ""

                for line in gfadata:
                    gfaall += line
                    line = line.strip('\n')
                    record = line.split('\t')
                    if record[0] == 'S':
                        seqlens.append(len(record[2]))

                newgfastore = GfaStore(run=runinstance, barcode = bar, readtype = ty)
                newgfastore.nreads = totfq
                newgfastore.gfaformat = gfaall  #   string  The whole GFA file
                newgfastore.save()

                #### SAVE 0 IF ASSSEMBLY FAILS

                newgfa = GfaSummary(run=runinstance, barcode = bar, readtype = ty)
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

    #try:
    if True:
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

        runidset = set()
        if inputtype == "flowcell":
            realflowcell = FlowCell.objects.get(pk=runid)
            flowcell_runs = FlowCellRun.objects.filter(flowcell=runid)
            for flowcell_run in flowcell_runs:
                runidset.add(flowcell_run.run_id)
                #print (flowcell_run.run_id)
                # we need to get the runids that make up this run
        else:
            runidset.add(runid)


        fastqs = FastqRead.objects.filter(run_id__id__in=runidset, id__gt=int(last_read))[:2000]

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
            runinstance = MinIONRun.objects.get(pk=runid)

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
                newpaf = PafStore(flowcell=realflowcell, read=readid, read_type=typeid)
                newpafstart = PafRoughCov(flowcell=realflowcell,read_type=typeid,barcode=fastqbarcode[record[0]],barcodegroup=fastqbarcodegroup[record[0]])
                newpafend = PafRoughCov(flowcell=realflowcell,read_type=typeid,barcode=fastqbarcode[record[0]],barcodegroup=fastqbarcodegroup[record[0]])
            elif inputtype == "run":
                newpaf = PafStore(run=run, read=readid, read_type=typeid)
            newpaf.reference = Reference

            logger.info("---> Before parsing paf record")
            logger.info(record)

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

        print('!!!!!!!It took {} to parse the paf.!!!!!!!!!'.format((donepafproc-doneminimaps)))

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
                                flowcell=realflowcell,
                                read_type=ty,
                                barcodegroup=bc,
                                reference=ref,
                                chromosome=ch,
                            )
                        summarycov.read_count += len(resultstore[ref][ch][bc][ty]['read'])
                        summarycov.cumu_length += resultstore[ref][ch][bc][ty]['length']
                        summarycov.save()


        jobdone = time.time()
        print('!!!!!!!It took {} to process the resultstore.!!!!!!!!!'.format((jobdone - donepafproc)))

    #except Exception as exception:
    #    print('An error occurred when running this task.')
    #    print(exception)
    #    JobMaster.objects.filter(pk=job_master_id).update(running=False)
    
    else:
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
            runinstance = MinIONRun.objects.get(pk=runid)
        elif inputtype == "flowcell":
            runinstance = FlowCell.objects.get(pk=runid)
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
    def __init__(self):
        self.tmpfile = tempfile.NamedTemporaryFile(suffix=".fa")
        self.krakenfile = tempfile.NamedTemporaryFile(suffix=".out")

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
        print('kraken-report --db /Volumes/SSD/kraken/minikraken_20141208 ' + self.krakenfile.name)
        p1 = subprocess.Popen('kraken-report --db /Volumes/SSD/kraken/minikraken_20141208 ' + self.krakenfile.name,
                              shell=True, stdout=subprocess.PIPE)
        (out, err) = p1.communicate()
        return out

    def process2(self):
        print('kraken-mpa-report --db /Volumes/SSD/kraken/minikraken_20141208 ' + self.krakenfile.name)
        p1 = subprocess.Popen('kraken-mpa-report --db /Volumes/SSD/kraken/minikraken_20141208 ' + self.krakenfile.name,shell=True, stdout=subprocess.PIPE)
        (out, err) = p1.communicate()
        return out

    def run(self):
        print('kraken --quick --db /Volumes/SSD/kraken/minikraken_20141208 --fasta-input --preload  ' + self.tmpfile.name)  # +'  #| kraken-translate --db /Volumes/SSD/kraken/minikraken_20141208/ $_'
        p1 = subprocess.Popen('kraken --quick --db /Volumes/SSD/kraken/minikraken_20141208 --fasta-input --preload  ' + self.tmpfile.name + '  ',shell=True, stdout=subprocess.PIPE)
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
                    runinstance = MinIONRun.objects.get(pk=runid)
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

    runs = MinIONRun.objects.all()

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
    MinIONRun.objects.filter(to_delete=True).delete()


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
