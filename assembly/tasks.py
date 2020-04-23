import gzip
import os
import subprocess
from datetime import timedelta, datetime
import tempfile
import numpy as np
import pandas as pd
from celery.task import task
from celery.utils.log import get_task_logger
from assembly.models import GfaStore, GfaSummary
from reads.utils import getn50
from reads.models import JobMaster, FastqRead, Flowcell, Run
logger = get_task_logger(__name__)


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


def callfetchreads_assem(runs,chunk_size,last_read):
    fastasm=list()
    #fastqs_list = list()
    fastq_df_barcode = pd.DataFrame()
    while True:
        #reads, last_read, read_count, fastasmchunk, fastqs_list_chunk = fetch_reads_alignment(runs, chunk_size, last_read)
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