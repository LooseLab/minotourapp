import pandas as pd
from celery.utils.log import get_task_logger

from reads.models import (
    FastqRead,
    FlowcellStatisticBarcode)

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


def create_flowcell_model(row,flowcell):
    """
    Create and return a series of model objects for bulk creation
    :param row: The row of the dataframe
    :type row: pandas.core.series.Series
    :return: A model objects to save into the dataframe
    """
    return FlowcellStatisticBarcode(
        flowcell=flowcell,
        barcode_name=row["barcode_name"],
        rejection_status=row["rejection_status"],
        read_type_name=row["read_type_name"],
        status=row["status"],
        total_length=row["total_length"],
        quality_sum=row["quality_sum"],
        read_count=row["read_count"],
        channel_count=row["channel_count"],
        iteration_count=row["iteration_count"],
        min_length=row["min_length"],
        max_length=row["max_length"],
        channel_presence=row["channel_presence"],
        sample_time=row["sample_time"]
    )