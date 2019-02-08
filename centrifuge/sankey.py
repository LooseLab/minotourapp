"""
sankey.py
"""
from centrifuge.models import CentrifugeOutput, SankeyLink
from jobs.models import JobMaster
import pandas as pd
from celery.utils.log import get_task_logger
import numpy as np

logger = get_task_logger(__name__)


def create_sankeylink_models(row):
    """
    apply to the dataframe and create sankeylinks objects for each row, return into a series
    :param row: The data frame row
    :return: The list of created objects
    """
    return SankeyLink(source=row["source"],
                      target=row["target"],
                      tax_id=row["tax_id"],
                      task=row["job_master"],
                      target_tax_level=row["target_tax_level"],
                      value=row["value"],
                      barcode_name=row["barcode_name"],
                      path=row["path"])


def calculate_sankey_values(lineages_df, flowcell, tax_rank_filter, metagenomics_task, barcode):
    """
    Calculate the values for the sankey diagram and save them into the database.
    :param lineages_df: A dataframe of lineages created from the current centrifuge output
    :param flowcell: The flowcell for these reads
    :param tax_rank_filter: The taxonomic ranks in descending order
    :param metagenomics_task: The task object relating to this metagenomics analysis
    :param barcode: The barcode that tjis sankey linkes object is running on
    :return:
    """
    # ##############  Storing snakey links and nodes ################

    logger.info("Flowcell id: {} - Starting Sankey Calculations".format(flowcell.id))

    # set index to tax_id
    lineages_df.set_index("tax_id", inplace=True)

    # The columns we actually need
    columns = ["num_matches", "sum_unique", "barcode_name", "superkingdom", "phylum", "class", "order", "family",
               "genus", "species"]

    # Subset the dataframe into those columns
    sankey_lineages_df = lineages_df[columns]

    # Backfill the dataframe nas with the lowest common entry, if family is missing, and genus
    #  isn't, family NaN becomes the genus entry
    sankey_lineages_df = sankey_lineages_df.fillna(axis=1, method="bfill")

    sankey_lineages_df = sankey_lineages_df.sort_values("num_matches", ascending=False)

    # Rank lineages to flow through
    sankey_lineages_df["path"] = np.arange(0, sankey_lineages_df.shape[0], 1)

    # initialise the series that will become our sankey pandas DataFrame
    source = pd.Series()
    target = pd.Series()
    value = pd.Series()
    tax_id = pd.Series()
    barcode = pd.Series()
    target_tax_level = pd.Series()
    path = pd.Series()
    # Iterate across the dataframe taxonnomic rank columns, to create the source, target, num reads series
    for index, rank in enumerate(tax_rank_filter):
        # The taxonomic rank of the source of this link
        source_tax_rank = rank
        # The next index up, so for kingdom it would be phylum
        shifted_index = index + 1
        if shifted_index < len(tax_rank_filter):
            # The target taxonomic rank
            target_tax_rank = tax_rank_filter[shifted_index]
            # Create a source series, append it to the bottom of the existing source series
            source = source.append(sankey_lineages_df[source_tax_rank])
            # Create a target series, append it to the bottom of the existing target series
            target = target.append(sankey_lineages_df[target_tax_rank])
            # Create a values series, append it to the bottom of the existing value series,
            # contains the number of matches
            value = value.append(sankey_lineages_df["num_matches"])
            # Create a tax_id series, append it to the bottom of the existing tax_id series
            tax_id = tax_id.append(pd.Series(sankey_lineages_df.index.values, index=sankey_lineages_df.index))
            # Create a barcode series, append it to the bottom of the existing barcode series
            barcode = barcode.append(sankey_lineages_df["barcode_name"])
            # Create a target tax level series, by broadcasting the value down the dataframe,
            # add it to the bottom of the existing series
            # Contains the taxonomic rank of the links target, useful for ordering later
            sankey_lineages_df["target_tax_level"] = target_tax_rank
            target_tax_level = target_tax_level.append(sankey_lineages_df["target_tax_level"])
            path = path.append(sankey_lineages_df["path"])

    # Create a DataFrame of links
    source_target_df = pd.concat([source, target, value, tax_id, target_tax_level, barcode, path], axis=1)

    # Rename the columns
    source_target_df.columns = ["source", "target", "value", "tax_id", "target_tax_level", "barcode_name", "path"]

    logger.info(source_target_df)

    # Create a series of Job_master objects through broadcasting
    source_target_df["job_master"] = metagenomics_task

    logger.info("Flowcell id: {} - Barcodes are:".format(flowcell.id))
    logger.info(source_target_df["barcode_name"].unique())

    # TODO DO sorting magic before the database deposition
    # Create a dataframe of the previous results for this metagenomics task

    previous = SankeyLink.objects.filter(task=metagenomics_task, barcode_name=barcode)
    if previous:
        previous.delete()
    # If there no results, i.e this is the first iteration

    to_create_sankey_df = source_target_df
    # if there are results these need updating where we have the same link

    # ## ### Create the SankeyLinks object #### # One for each link ## ### #
    to_create_sankey_df.reset_index(inplace=True)

    # Drop duplicates as we only need to create one
    to_create_sankey_df = to_create_sankey_df.drop_duplicates(subset=["tax_id", "target_tax_level"])

    # If there are values that need creating
    if not to_create_sankey_df.empty:
        # Perform the apply to create the series to create
        sankey_models_for_bulk_insert = to_create_sankey_df.apply(create_sankeylink_models, axis=1)

        # Bulk create the objects
        SankeyLink.objects.bulk_create(list(sankey_models_for_bulk_insert.values))

    logger.info("Flowcell id: {} - Finished updating sankey links".format(flowcell.id))


def calculate_sankey(flowcell_job_id):
    """
    Run the calculate sankey task
    :param flowcell_job_id: The id (pk) of the jobmaster object
    :return:
    """
    task = JobMaster.objects.get(pk=flowcell_job_id)
    flowcell = task.flowcell
    task.running = True
    task.save()
    metagenomics_task = JobMaster.objects.get(flowcell=flowcell, job_type__name="Metagenomics")
    tax_rank_filter = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    barcodes = CentrifugeOutput.objects.filter(task=metagenomics_task).values_list("barcode_name", flat=True).distinct()

    for barcode in barcodes:

        df = pd.DataFrame(list(CentrifugeOutput.objects.filter(
            task=metagenomics_task, barcode_name=barcode).order_by("-num_matches").values()[:30]))
        df.rename(columns={"classy": "class"}, inplace=True)
        logger.info(df.head())

        calculate_sankey_values(df, flowcell, tax_rank_filter, metagenomics_task, barcode)
    task.running = False
    task.complete = True
    task.save()
