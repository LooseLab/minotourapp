"""contains utility functions for the  centrifuge.py file """

import numpy as np
import pandas as pd

from metagenomics.models import CentrifugeOutput, DonutData, Metadata


def convert_species_to_subspecies(name):
    """
    Centrifuge doesn't output subspecies as a rank. Checks if subsp. is in a species name, and creates the
    corresponding rank and subspecies name
    Parameters
    ----------
    name - (type string) The name of the strain

    Returns - (type String) The subspecies as a string
    -------

    """
    # if the name is a string, not a float
    if type(name) is str:
        # Fairly obvious
        if "subsp." in name:
            basestring = name.split("subsp.", 1)[0]
            string = name.split("subsp.", 1)[1]
            string = string.split(" ")[1]
            subspecies = basestring + " subsp. " + string
            return subspecies
    else:
        pass


def delete_series(series, df):
    """

    Parameters
    ----------
    series : (type list of Series) - series to be deleted
    df : (type a pandas DataFrame) - DataFrame to delete them from

    :return : a DataFrame with the series in the list deleted
    -------

    """
    for s in series:
        if type(df[s]).__name__ == "Series":
            df = df.drop(s, 1)
    return df


def divd(row, cl_bar):
    """
    Calculate_proportion_of_classified
    :param row: The df row
    :param cl_bar: The number of reads classified in a barcode, in dict form keyed to the barcode_name
    :return:
    """
    return round(
        (row["num_matches"] / cl_bar[row["barcode_name"]]) * 100, 4
    )


def create_centrifuge_models(row, classified_per_barcode):
    """
    Append a CentOutput object to a list for each row in the metagenomics output
    :param row: the row from the data frame
    :param classified_per_barcode: The number of reads classified for each barcode as a dictionary
    :param first: Boolean saying if this is the first iteration, if it isn't a pandas merge adds a suffix onto the code
    :return: The list of newly created objects

    """
    # If there is a NaN in the proportion of classified, redo the calculation
    if (
        row["proportion_of_classified"] == "Unclassified"
        or row["proportion_of_classified"] == np.NaN
    ):
        row["proportion_of_classified"] = round(
            row["num_matches"] / classified_per_barcode[row["barcode_name"]], 3
        )

    return CentrifugeOutput(
        name=row["name"],
        tax_id=row["tax_id"],
        task=row["task"],
        num_matches=row["num_matches"],
        sum_unique=row["sum_unique"],
        barcode_name=row["barcode_name"],
        proportion_of_classified=row["proportion_of_classified"],
        superkingdom=row["superkingdom"],
        phylum=row["phylum"],
        classy=row["class"],
        order=row["order"],
        family=row["family"],
        genus=row["genus"],
        species=row["species"],
        latest=row["latest"],
    )


def create_donut_data_models(row, task):
    """
    Create and return a series of model objects for bulk creation
    :param row: The row of the dataframe
    :type row: pandas.core.series.Series
    :param task: The task objects for this analysis
    :return: A model objects to save into the dataframe
    """

    return DonutData(
        task=task,
        num_matches=row["num_matches"],
        sum_unique=row["sum_unique"],
        name=row["name"],
        tax_rank=row["tax_rank"],
        barcode_name=row["barcode_name"],
        latest=row["latest"],
    )


def falls_in_region(row, map_df):
    """
    Does this reads mapping fall within or encompass any of the regions we have defined?
    :param row: The gff3 dataframe row
    :param map_df_by_species: The results of the mapping step in a pandas dataframe
    :return:
    """
    # TODO include which region it maps to
    # If the start is higher than the region start
    map_df_by_species = map_df[map_df["name"] == row["species"].replace(" ", "_")]
    start_low = map_df_by_species["target_start"] > row["start"]
    # If the start is lower than the region end
    start_high = map_df_by_species["target_start"] < row["end"]
    # If the end is higher than the region end
    end_low = map_df_by_species["target_end"] > row["start"]
    # If the end is lower than the region end
    end_high = map_df_by_species["target_end"] < row["end"]
    # If the read completely spans the gene/region
    read_is_wider_than_region = ((map_df_by_species["target_start"] < row["start"]) & (map_df_by_species["target_end"] > row["end"]))
    # Concat the series together, create "truth" dataframe
    map_bool_df = pd.concat(
        [start_low, start_high, end_low, end_high, read_is_wider_than_region], axis=1, ignore_index=True
    )
    # Keep where the start is in the region,
    map_bool_df["keep"] = np.where(
        (map_bool_df[0] & map_bool_df[1]) | (map_bool_df[2] & map_bool_df[3] | map_bool_df[4]), 1, 0
    )
    return map_bool_df["keep"].any()


def calculate_proportion_for_table(row, task):
    """

    Parameters
    ----------
    row: pandas.core.series.Series
        Dataframe row
    task: reads.models.JobMaster
        Task orm from the centrifuge job
    Returns
    -------
    list
    """
    species_as_percent_of_total = round(row["num_matches"]/Metadata.objects.get(task=task).classified*100, 2)
    mapped_as_percent_of_classed = round(row["num_mapped"]/row["num_matches"]*100, 2)
    mapped_to_target_area_percent = round(row["red_reads"]/row["num_mapped"]*100, 2)
    return [species_as_percent_of_total, mapped_as_percent_of_classed, mapped_to_target_area_percent]


def get_metagenomics_data(task, flowcell):
    """
    Get metagenomics output for given task
    Parameters
    ----------
    task: reads.models.JobMaster
        Metagenomics JobMaster
    flowcell: minknow_data.models.Flowcell
        Flowcell that run was on

    Returns
    -------
    pandas.DataFrame
    """
    centrifuge_df = pd.DataFrame(
        list(
            CentrifugeOutput.objects.filter(task=task)
            .exclude(barcode_name="No")
            .values()
        )
    )
    centrifuge_df.rename(columns={"classy": "class"}, inplace=True)
    column_order = [
        "barcode_name",
        "tax_id",
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "num_matches",
        "proportion_of_classified",
    ]
    return centrifuge_df[column_order]