import numpy as np
import pandas as pd
from celery import task
from celery.utils.log import get_task_logger

from django.db.models import Avg, Max, Min, Sum, Q

from jobs.models import JobMaster
from reads.models import (
    FastqRead,
    HistogramSummary,
    Flowcell,
    FlowcellSummaryBarcode,
    FlowcellStatisticBarcode,
    FlowcellHistogramSummary,
    FlowcellChannelSummary,
)
from reads.services import (
    save_flowcell_summary_barcode,
    save_flowcell_statistic_barcode,
    save_flowcell_histogram_summary,
    save_flowcell_channel_summary,
)
from readuntil.models import RejectedFastqRead, AcceptedFastqRead
from reads.models import Barcode


logger = get_task_logger(__name__)


def callfetchreads(runs, chunk_size, last_read):
    fastq_df_barcode = pd.DataFrame()
    while True:
        reads, last_read, read_count = fetchreads(runs, chunk_size, last_read)
        fastq_df_barcode = fastq_df_barcode.append(reads)
        if len(fastq_df_barcode) >= chunk_size or len(reads) == 0:
            break
    read_count = len(fastq_df_barcode)
    return fastq_df_barcode, last_read, read_count


def fetchreads(runs, chunk_size, last_read):
    countsdict = dict()
    for run in runs:
        fastqs = (
            FastqRead.objects.values_list("id")
            .filter(run=run, id__gt=int(last_read))
            .first()
        )
        if fastqs:
            countsdict[fastqs[0]] = run
    count = 1
    fastq_df_barcode = pd.DataFrame()
    if len(countsdict) > 1:
        for entry in sorted(countsdict):
            if count < len(countsdict):
                fastqs = FastqRead.objects.filter(
                    run=countsdict[entry],
                    id__gt=int(last_read),
                    id__lt=int(sorted(countsdict)[count]),
                )[:chunk_size]
                fastq_df_barcode = fastq_df_barcode.append(
                    pd.DataFrame.from_records(
                        fastqs.values(
                            "id",
                            "start_time",
                            "barcode__name",
                            "rejected_barcode__name",
                            "type__name",
                            "is_pass",
                            "sequence_length",
                            "quality_average",
                            "channel",
                        )
                    )
                )

                last_read = fastq_df_barcode.id.max()

                if len(fastq_df_barcode) >= chunk_size:
                    break
            count += 1
    elif len(countsdict) == 1:
        """This is risky and it breaks the logic - we end up grabbing reads"""
        mykey = list(countsdict.keys())[0]
        fastqs = FastqRead.objects.filter(
            run=countsdict[mykey], id__gt=int(last_read),
        )[:chunk_size]
        fastq_df_barcode = fastq_df_barcode.append(
            pd.DataFrame.from_records(
                fastqs.values(
                    "id",
                    "start_time",
                    "barcode__name",
                    "rejected_barcode__name",
                    "type__name",
                    "is_pass",
                    "sequence_length",
                    "quality_average",
                    "channel",
                )
            )
        )
        last_read = fastq_df_barcode.id.max()

    read_count = len(fastq_df_barcode)
    return fastq_df_barcode, last_read, read_count

#
# def compare_read_ids(flowcell, runs):
#     """
#     Compare read ids in the rejected fastqread table and attach a barcode to the rejected read ID table
#     :param flowcell: The flowcell the read ids are attached to
#     :type flowcell: reads.models.Flowcell
#     :param runs: A queryset of the runs so we can create the rejected/accepted barcodes
#     :type runs: django.db.models.query.QuerySet
#     :return:
#     """
#
#     # A list of tuples, with the model for the read_ids to be looked up
#     # and the name of the barcode to add to the FastQRead
#     tupley_listy = [(RejectedFastqRead, "U"), (AcceptedFastqRead, "S")]
#
#     # loop through the runs
#     for run in runs:
#         for read_type in tupley_listy:
#             # get the rejected / accepted barcode
#             rejected_barcode, created = Barcode.objects.get_or_create(
#                 run=run, name=read_type[1]
#             )
#             # Get the readIds from the rejected / accepted tables
#             rejected_ids = (
#                 read_type[0]
#                 .objects.filter(flowcell=flowcell)
#                 .values_list("read_id", flat=True)
#             )
#             # update the matching fastqread objects with whichever type of Read the for loop
#             # is on (Accepted or Rejected)
#             FastqRead.objects.exclude(rejected_barcode=rejected_barcode).filter(
#                 run=run, read_id__in=rejected_ids
#             ).update(barcode=rejected_barcode)


@task()
def chancalc(flowcell_id, job_master_id, last_read):
    """
    Basecalled data metrics calculations
    :param flowcell_id: Id of the flowcell for these reads
    :type flowcell_id: int
    :param job_master_id: The id of the chancalc JobMaster database entry
    :type job_master_id: int
    :param last_read: The last read we have looked at in previous chancalc iterations. 0 if this is the first time.
    :type last_read: int
    :return: None
    """

    job_master = JobMaster.objects.get(pk=job_master_id)

    flowcell = Flowcell.objects.get(pk=flowcell_id)

    runs = flowcell.runs.all()

    compare_read_ids(flowcell, runs)

    chunk_size = 50000
    # Get the fastq data
    fastq_df_barcode, last_read, _ = callfetchreads(runs, chunk_size, last_read)

    fastqlen = fastq_df_barcode.shape[0]

    if fastqlen > 0:

        new_last_read = fastq_df_barcode.iloc[-1]["id"]
        #
        # if read_count and last read are 0, then delete any previous summaries summaries
        #
        if job_master.read_count == 0 and job_master.last_read == 0:

            print("Flowcell id: {} - Deleting summaries".format(flowcell.id))

            FlowcellSummaryBarcode.objects.filter(flowcell=flowcell).delete()
            FlowcellStatisticBarcode.objects.filter(flowcell=flowcell).delete()
            FlowcellHistogramSummary.objects.filter(flowcell=flowcell).delete()
            FlowcellChannelSummary.objects.filter(flowcell=flowcell).delete()

        fastq_df_barcode["status"] = np.where(
            fastq_df_barcode["is_pass"] == False, "Fail", "Pass"
        )

        fastq_df_barcode["start_time_truncate"] = np.array(
            fastq_df_barcode["start_time"], dtype="datetime64[m]"
        )

        print(
            "Flowcell id: {} - The new last read is {}".format(
                flowcell.id, new_last_read
            )
        )

        fastq_df_allreads = fastq_df_barcode.copy()
        fastq_df_allreads["barcode__name"] = "All reads"

        fastq_df = fastq_df_barcode.append(fastq_df_allreads)
        print(f"fastq_df_barcode is {fastq_df_barcode.head()}")
        print(f"fastq_df is {fastq_df.head()}")

        #
        # Calculates statistics for flowcellSummaryBarcode
        #
        fastq_df_result = fastq_df.groupby(
            [
                "barcode__name",
                "rejected_barcode__name",
                "type__name",
                "is_pass",
            ]
        ).agg(
            {
                "sequence_length": ["min", "max", "sum", "count"],
                "quality_average": ["sum"],
                "channel": ["unique"],
            }
        )
        logger.info(fastq_df_result.head())
        logger.info(f"keys is {fastq_df_result.keys()}")
        fastq_df_result.reset_index().apply(
            lambda row: save_flowcell_summary_barcode(flowcell_id, row), axis=1
        )

        # fastq_df['start_time']=fastq_df['start_time'].values.astype('<M8[m]')
        #
        # Calculates statistics for FlowcellStatisticsBarcode
        #
        fastq_df_result = fastq_df.groupby(
            [
                "start_time_truncate",
                "barcode__name",
                "type__name",
                "is_pass",
                "rejected_barcode__name",
            ]
        ).agg(
            {
                "sequence_length": ["min", "max", "sum", "count"],
                "quality_average": ["sum"],
                "channel": ["unique"],
            }
        )

        fastq_df_result.reset_index().apply(
            lambda row: save_flowcell_statistic_barcode(flowcell_id, row),
            axis=1,
        )

        #
        # Calculates statistics for HistogramSummary
        #
        fastq_df["bin_index"] = (
            fastq_df["sequence_length"]
            - fastq_df["sequence_length"] % HistogramSummary.BIN_WIDTH
        ) / HistogramSummary.BIN_WIDTH

        fastq_df_result = fastq_df.groupby(
            [
                "barcode__name",
                "type__name",
                "is_pass",
                "bin_index",
                "rejected_barcode__name",
            ]
        ).agg({"sequence_length": ["sum", "count"]})

        fastq_df_result.reset_index().apply(
            lambda row: save_flowcell_histogram_summary(flowcell_id, row),
            axis=1,
        )

        #
        # Calculates statistics for ChannelSummary
        #
        fastq_df_result = fastq_df.groupby(["channel"]).agg(
            {"sequence_length": ["sum", "count"]}
        )

        fastq_df_result.reset_index().apply(
            lambda row: save_flowcell_channel_summary(flowcell_id, row), axis=1
        )

        # last_read = fastq_df_barcode['id'].max()
        last_read = new_last_read

    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = False
    job_master.last_read = last_read
    job_master.read_count = job_master.read_count + fastqlen
    job_master.save()
