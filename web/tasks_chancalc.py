import numpy as np
import pandas as pd
from celery import task
from celery.utils.log import get_task_logger

from jobs.models import JobMaster
from reads.models import FastqRead, HistogramSummary, Flowcell, FlowcellSummaryBarcode, FlowcellStatisticBarcode, \
    FlowcellHistogramSummary, FlowcellChannelSummary
from reads.services import save_flowcell_summary_barcode, save_flowcell_statistic_barcode, \
    save_flowcell_histogram_summary, save_flowcell_channel_summary

logger = get_task_logger(__name__)


@task()
def chancalc(flowcell_id, job_master_id, last_read):

    job_master = JobMaster.objects.get(pk=job_master_id)
    # job_master.running = True
    # job_master.save()

    flowcell = Flowcell.objects.get(pk=flowcell_id)

    ## This is a slow query and needs optimisation...

    #fastqs = FastqRead.objects.filter(run__flowcell=flowcell).filter(id__gt=int(last_read)).order_by('id')[:50000]

    runs = flowcell.runs.all()

    chunk_size = 20000

    holderformaxreadseen = 0

    for run in runs:
        # This query doesn't work because it assumes that reads are entered into the database in run order
        # And they are not. So we need a way to work around this.
        #
        '''
        fastqs = FastqRead.objects.filter(run=run, id__gt=int(task.last_read)
                                          ).order_by('id')[:chunk_size]

        if fastqs.count() > 0:
            # we have data - so go off and use it.
            proceed = True
            break
        '''
        if holderformaxreadseen > 0:
            fastqs = FastqRead.objects.filter(run=run, id__gt=int(job_master.last_read),
                                              id__lt=int(holderformaxreadseen)).order_by('id')
        else:
            fastqs = FastqRead.objects.filter(run=run, id__gt=int(job_master.last_read)).order_by('id')[:chunk_size]

        if 'fastq_store' in locals():
            fastq_store = fastq_store.union(fastqs)
        else:
            fastq_store = fastqs

        if holderformaxreadseen == 0 and fastqs.count() > 0:
            holderformaxreadseen = fastqs[fastqs.count() - 1].id

    # May not always create fastqstore need to check.
    if runs.count() > 1:
        fastqs = fastq_store.order_by('id')[:chunk_size]

    fastq_df_barcode = pd.DataFrame.from_records(
        fastqs.values('id', 'start_time', 'barcode__name', 'type__name', 'is_pass', 'sequence_length',
                      'quality_average', 'channel'))

    fastqlen = len(fastq_df_barcode)

    if fastqlen > 0:
        new_last_read = fastq_df_barcode.iloc[-1]['id']
        #
        # if read_count = 0, then delete summaries
        #
        if job_master.read_count == 0 and job_master.last_read == 0:

            logger.info('Flowcell id: {} - Deleting summaries'.format(flowcell.id))

            FlowcellSummaryBarcode.objects.filter(flowcell=flowcell).delete()
            FlowcellStatisticBarcode.objects.filter(flowcell=flowcell).delete()
            FlowcellHistogramSummary.objects.filter(flowcell=flowcell).delete()
            FlowcellChannelSummary.objects.filter(flowcell=flowcell).delete()

        #last_fastq = fastqs[len(fastqs) - 1]

        #new_last_read = last_fastq.id


        #fastq_df_barcode = pd.DataFrame.from_records(fastqs.values('id', 'start_time', 'barcode__name', 'type__name', 'is_pass', 'sequence_length', 'quality_average', 'channel'))

        fastq_df_barcode['status'] = np.where(fastq_df_barcode['is_pass'] == False, 'Fail', 'Pass')
        fastq_df_barcode['start_time_truncate'] = np.array(fastq_df_barcode['start_time'], dtype='datetime64[m]')

        logger.info('Flowcell id: {} - The new last read is {}'.format(flowcell.id, new_last_read))

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

        # last_read = fastq_df_barcode['id'].max()
        last_read = new_last_read


    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = False
    job_master.last_read = last_read
    job_master.read_count = job_master.read_count + fastqlen
    job_master.save()