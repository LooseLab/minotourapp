"""contains utility functions for the  centrifuge.py file """
from collections import defaultdict

import numpy as np
import pandas as pd
from ete3 import NCBITaxa

from metagenomics.models import MappingResult, CentrifugeOutput, DonutData


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


def create_mapping_results_objects(barcodes, task, temp_targets):
    """
    Create the mapping results objects for the
    Parameters
    ----------
    barcodes: list
        The list of barcodes that we have to create the mapping results objects for
    task: reads.models.JobMaster
        The django ORM task of this metagenomics task
    temp_targets: list
        The list of targets we have

    Returns
    -------

    """
    ncbi = NCBITaxa()
    already_created_barcodes = list(
        MappingResult.objects.filter(task=task)
        .values_list("barcode_name", flat=True)
        .distinct()
    )
    # for each barcode in the target_df
    for barcode in barcodes:
        # Check if barcode has mapping result objects created for it already
        if barcode not in already_created_barcodes:
            # for each target, create a Mapping result
            for target in temp_targets:
                # get the tax_id for this target
                tax_id = ncbi.get_name_translator([target])[target][0]
                mr = MappingResult(
                    task=task, barcode_name=barcode, species=target, tax_id=tax_id,
                )
                mr.save()


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


def calculate_lineages_df(ncbi, df, tax_rank_filter, flowcell):
    """
    Insert any new lineages of species that are in our results into the database, and fetch all that we need for these
    results
    :param ncbi: The instance of the NCBITaxa class
    :param df: The data frame of our results
    :param tax_rank_filter: The order of our taxonomic ranks
    :param flowcell: The flowcell object
    :return lineages_Df: a DataFrame of all the lineages
    """
    # create new defaultDict which creates a default dict for missing values
    lineage_dict = defaultdict(lambda: defaultdict())
    # Get the current taxIDs in the Dataframe
    taxid_list = df["tax_id"].values
    # Check the taxIDs we have in the dataframe, to see if they have matching lineages already.
    # If not the set subtraction means that unique, lineage-less taxIDs end up in the not_prev_lineaged dict
    # Get the lineages for these taxIDs in the new metagenomics output,
    # returns a dict keyed by taxID,
    # value is listof taxIDS in lineage
    lineages_taxidlist_dict = ncbi.get_lineage_translator(taxid_list)
    # loop over the new lineages dictionary, the values of which is a list of taxIDs from root to the organism,
    # key is a int taxID
    # TODO rewrite into pandas mapping function, low priority, or magic datafame magic
    for key, value in lineages_taxidlist_dict.items():
        # get the ranks of a list of taxID. Returns a dict with keys of the taxID in the list of taxIDS,
        # value is the rank, i.e root, kingdom etc.
        # TODO vectorise this sucka
        lineage_ranked = ncbi.get_rank(value)
        # a dict where key is taxID, value is the name i.e {562: Escherichia coli}
        taxid_species_lookup_dict = ncbi.get_taxid_translator(value)
        # Populate the defaultDict, by looping the lineage ranked dict and adding it where they fit the
        # condition, the if key==taxID is if they are subspecies or strain
        lineage_dict[key] = {
            rank: taxid_species_lookup_dict[tax_id]
            for tax_id, rank in lineage_ranked.items()
            if rank in tax_rank_filter or key == tax_id
        }
    # create dataframe, the .T means that the dataframe is transposed, so dict key is the series title
    lineages_df = pd.DataFrame.from_records(lineage_dict).T
    # rename the no rank column to strain
    lineages_df = lineages_df.rename(columns={"no rank": "strain"})
    # if these series don't exist add them and populate with numpy NaN so there's no key errors
    for series in ["subspecies", "strain", "subStrainSpecies"]:
        if series not in lineages_df.keys():
            lineages_df[series] = np.NaN
    print("Flowcell id: {} - Determining Subspecies".format(flowcell.id))
    # create new additional subspecies column, from strains column which has subspecies above it in taxa level
    lineages_df["subStrainSpecies"] = lineages_df["strain"].map(
        convert_species_to_subspecies
    )
    # merge new subspecies column with existing column
    lineages_df["subspecies"].fillna(
        lineages_df["subStrainSpecies"], inplace=True
    )
    unclassified_row = {key: np.NaN for key in lineages_df.keys()}
    row = pd.DataFrame(unclassified_row, index=[0])
    lineages_df = lineages_df.append(row, sort=True)
    # delete the new subspecies column
    delete_series(["subStrainSpecies"], lineages_df)
    return lineages_df
