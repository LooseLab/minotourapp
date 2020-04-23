import numpy as np
import pandas as pd
from celery import task
from celery.utils.log import get_task_logger


from reads.models import (
    FastqRead,
    HistogramSummary,
    Flowcell,
    FlowcellSummaryBarcode,
    FlowcellStatisticBarcode,
    FlowcellHistogramSummary,
    FlowcellChannelSummary,
    JobMaster)
from reads.services import (
    save_flowcell_summary_barcode,
    save_flowcell_statistic_barcode,
    save_flowcell_histogram_summary,
    save_flowcell_channel_summary,
)


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

    #last_iteration = job_master.iteration_count

    flowcell = Flowcell.objects.get(pk=flowcell_id)

    runs = flowcell.runs.all()

    chunk_size = 500000
    # Get the fastq data
    fastq_df_barcode, last_read, _ = callfetchreads(runs, chunk_size, last_read)

    fastqlen = fastq_df_barcode.shape[0]

    #
    # if read_count and last read are 0, then delete any previous summaries summaries
    #
    if job_master.read_count == 0 and job_master.last_read == 0:
        print("Flowcell id: {} - Deleting summaries".format(flowcell.id))

        FlowcellSummaryBarcode.objects.filter(flowcell=flowcell).delete()
        FlowcellStatisticBarcode.objects.filter(flowcell=flowcell).delete()
        FlowcellHistogramSummary.objects.filter(flowcell=flowcell).delete()
        FlowcellChannelSummary.objects.filter(flowcell=flowcell).delete()

    if fastqlen > 0:

        new_last_read = fastq_df_barcode.iloc[-1]["id"]


        fastq_df_barcode["status"] = np.where(
            fastq_df_barcode["is_pass"] == False, "Fail", "Pass"
        )

        fastq_df_barcode["start_time_truncate"] = np.array(
            fastq_df_barcode["start_time"], dtype="datetime64[m]"
        )

        #print(
        #    "Flowcell id: {} - The new last read is {}".format(
        #        flowcell.id, new_last_read
        #    )
        #)

        fastq_df_allreads = fastq_df_barcode.copy()
        fastq_df_allreads["barcode__name"] = "All reads"

        fastq_df = fastq_df_barcode.append(fastq_df_allreads)

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
        logger.debug(fastq_df_result.head())
        logger.debug(f"keys is {fastq_df_result.keys()}")
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

        ### This is really slow - and I don't know why!

        #temp_hold = fastq_df_result.reset_index(inplace=True)
        # fastq_df_result.columns = fastq_df_result.columns.get_level_values(0)
        # #temp_hold['start_time'] = pd.to_datetime(temp_hold["start_time_truncate"])
        #BarcodeObjects = FlowcellStatisticBarcode.objects.filter(flowcell=flowcell)
        #tempframe = pd.DataFrame.from_records(BarcodeObjects.values())

        ### Manipulate the data to get the values we need in a specific data frame.

        barcoderesults = fastq_df_result.reset_index()
        barcoderesults['start_time'] = pd.to_datetime(barcoderesults["start_time_truncate"],utc=True)
        barcoderesults["sequence_length_sum"] = barcoderesults["sequence_length"]["sum"]
        barcoderesults["sequence_length_min"] = barcoderesults["sequence_length"]["min"]
        barcoderesults["sequence_length_max"] = barcoderesults["sequence_length"]["max"]
        barcoderesults["quality_average_sum"] = barcoderesults["quality_average"]["sum"]
        barcoderesults["read_count"] = barcoderesults["sequence_length"]["count"]
        barcoderesults["channels"] = barcoderesults["channel"]["unique"]
        barcoderesults["channel_count"] = len(barcoderesults["channels"])
        channel_presence = np.zeros(3000)
        idx = np.array(barcoderesults["channels"])-1
        channel_presence[idx[0]]=1
        barcoderesults["channel_presence"]="".join(np.char.mod('%d', channel_presence))
        testbarcoderesults = barcoderesults
        testbarcoderesults = testbarcoderesults.drop(['sequence_length','quality_average','channel','start_time_truncate'],axis=1)
        testbarcoderesults.columns = testbarcoderesults.columns.get_level_values(0)
        testbarcoderesults = testbarcoderesults.rename(columns={"start_time": "sample_time", "barcode__name": "barcode_name", "type__name": "read_type_name", 'is_pass': 'status', 'rejected_barcode__name': 'rejection_status','sequence_length_sum': 'total_length', 'sequence_length_min': 'min_length', "sequence_length_max": 'max_length', 'quality_average_sum':'quality_sum'})
        testbarcoderesults["status"]=testbarcoderesults["status"].astype(str)

        #But we want to get the current set of results so:
        tempframe = pd.DataFrame.from_records(FlowcellStatisticBarcode.objects.filter(flowcell=flowcell,iteration_count=job_master.iteration_count).all().values())
        ###tempframe is the existing results
        if tempframe.empty:
            updateset=testbarcoderesults
        else:
            #This query merges existing records with new records
            mergeset = testbarcoderesults.merge(tempframe, how='outer', on=["barcode_name","read_type_name","sample_time","rejection_status","status"])
            #Now we merge the various results to get the data we need.
            mergeset["channel_presence_x"] = mergeset["channel_presence_x"].fillna('0'*3000)
            mergeset["channel_presence_y"] = mergeset["channel_presence_y"].fillna('0'*3000)
            mergeset = mergeset.fillna(0)
            mergeset["total_length"]=mergeset["total_length_x"]+mergeset["total_length_y"]
            mergeset["read_count"]=mergeset["read_count_x"]+mergeset["read_count_y"]
            mergeset["quality_sum"]=mergeset["quality_sum_x"]+mergeset["quality_sum_y"]

            #todo ##minimum value if not zero
            mergeset["min_length"] = mergeset[["min_length_x", "min_length_y"]].min(axis=1)
            mergeset["max_length"] = mergeset[["max_length_x", "max_length_y"]].max(axis=1)
            mergeset["channel_presence_x"]=mergeset["channel_presence_x"].apply(np.fromiter, dtype="int")
            mergeset["channel_presence_y"]=mergeset["channel_presence_y"].apply(np.fromiter, dtype="int")
            mergeset["channel_presence"]=np.bitwise_or(mergeset["channel_presence_x"],mergeset["channel_presence_y"])
            mergeset["channel_presence"]=mergeset["channel_presence"].apply(lambda x: x.astype("str")).str.join("")
            #mergeset["channels"]=set(pd.concat([mergeset["channels_x"], mergeset["channels_y"]]))
            mergeset["channel_count"]=mergeset['channel_presence'].str.count(pat="1")
            updateset = mergeset.drop(['channels', 'total_length_x', 'min_length_x', 'max_length_x','quality_sum_x', 'read_count_x', 'channel_presence_x','total_length_y', 'read_count_y', 'max_length_y','min_length_y', 'quality_sum_y', 'channel_presence_y'],axis=1)

        ### Now we change the iteration count.

        updateset["iteration_count"]=job_master.iteration_count + 1

        ### Now we create a set of records to insert into the database.

        to_bulk_save_series = updateset.apply(
            create_flowcell_model, args=(flowcell, ), axis=1
        )

        # Bulk create all those objects in the series
        FlowcellStatisticBarcode.objects.bulk_create(
            to_bulk_save_series.values.tolist(), batch_size=1000
        )

        job_master.iteration_count+=1
        job_master.save()

        ### DELETE OLD FLOWCELL STATISTIC BARCODES
        FlowcellStatisticBarcode.objects.filter(flowcell=flowcell, iteration_count=job_master.iteration_count-1).delete()
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
