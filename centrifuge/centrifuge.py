"""Centrifuge.py
"""
import subprocess
from collections import defaultdict
from io import StringIO
import numpy as np
import pandas as pd
from celery.utils.log import get_task_logger
from django.db.models import ObjectDoesNotExist
from django.utils import timezone
from ete3 import NCBITaxa
from centrifuge.models import CentrifugeOutput, LineageValue, Metadata, SankeyLink, CentrifugeOutputBarcoded, \
    MappingResult, TargetMappedReadId, MappingTarget, MappingResultsBarcoded
from jobs.models import JobMaster
from minotourapp.utils import get_env_variable
from reads.models import FastqRead
from reference.models import ReferenceInfo
from django.conf import settings
from django.db.utils import IntegrityError

pd.options.mode.chained_assignment = None
logger = get_task_logger(__name__)


# TODO del unused df and things

def calculate_barcoded_values(barcode_group_df, barcode_df):
    """
    Return the main cent output data frame with the barcode only result concatenated on
    :param barcode_group_df: barcode data frame group data frame, one for each barcode
    :param barcode_df: the main results data frame
    :return: A new df with the result for barcode appended
    """

    barcode_group_df_gb = barcode_group_df.groupby("tax_id")
    # calculate total unique matches for this species in this barcode as a series
    sum_unique = barcode_group_df_gb["unique"].sum()
    # calculate the total matches fot this species in this barcode as a series
    num_matches = barcode_group_df_gb.size()

    barcode_group_df.drop(columns=["readID", "seqID", "numMatches", "unique"], inplace=True)

    barcode_group_df.drop_duplicates(subset=["barcode_name", "name"], keep="first", inplace=True)

    # combine the dataframes
    values_df = pd.concat([sum_unique, num_matches, barcode_group_df["name"]], axis=1)

    barcode_df = barcode_df.append(values_df, sort=True)

    return barcode_df


def create_centrifuge_models(row):
    """
    Append a CentOutput object to a list for each row in the centrifuge output
    :param row: the row from the data frame
    :return: The list of newly created objects
    """
    return CentrifugeOutput(name=row["name"],
                            tax_id=row["tax_id"],
                            task=row["task"])


def create_centrifuge_barcoded_models(row, job_master):
    """
    Create centrifuge_output_barcoded objects and return as an element in a new series
    :param job_master: The job_master object for this task run
    :param row: The dataframe row
    :return: The  newly created objects
    """
    try:
        centrifuge_output = CentrifugeOutput.objects.get(tax_id=row["tax_id"], task__id=job_master.id)
        return CentrifugeOutputBarcoded(num_matches=row["num_matches"],
                                        sum_unique=row["sum_unique"],
                                        output=centrifuge_output,
                                        tax_id=row["tax_id"],
                                        barcode_name=row["barcode_name"])

    except ObjectDoesNotExist:
        logger.info("<<<<<")
        logger.info("Cent Object doesn't exist")
        logger.info("No matching entry for tax_id {} and task_id {}".format(row["tax_id"], job_master.id))
        logger.info("<<<<<")


def update_centrifuge_barcoded_values(row, flowcell_job_id):
    """
    Update existing CentrifugeOutputBarcoded objects in the database
    :param flowcell_job_id: The id for the job_master for this task run
    :param row: The data frame row
    :return: The list of newly created objects
    """
    CentrifugeOutputBarcoded.objects.filter(tax_id=row["tax_id"],
                                            barcode_name=row["barcode_name"],
                                            output__task__id=flowcell_job_id
                                            ).update(
        num_matches=row["updated_num_matches"],
        sum_unique=row["updated_sum_unique"]
    )


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

        if "subsp." in name:
            basestring = name.split("subsp.", 1)[0]
            string = name.split("subsp.", 1)[1]
            string = string.split(" ")[1]
            subspecies = basestring + " subsp. " + string
            return subspecies
    else:
        pass


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
                      barcode_name=row["barcode_name"])


def update_sankeylink_values(row, flowcell_job_id):
    """
    Update the barcode sankey links rows in the database
    :param row: The data frame row
    :param flowcell_job_id: The Pk of the task ID
    :return: Nothing
    """
    SankeyLink.objects.filter(task__id=flowcell_job_id,
                              target_tax_level=row["target_tax_level"],
                              barcode_name=row["barcode_name"],
                              tax_id=row["tax_id"]).update(value=row["updated_value"])


def create_lineage_models(row):
    """
    Apply function to the dataframe to populate a list with a LineageValue for each row
    :param row: The row of the dataframe
    :return: The list of newly created objects
    """
    return LineageValue(superkingdom=row["superkingdom"], phylum=row["phylum"],
                        tax_id=row["tax_id"],
                        classy=row["class"],
                        order=row["order"],
                        family=row["family"],
                        genus=row["genus"],
                        species=row["species"],
                        subspecies=row["subspecies"],
                        strain=row["strain"]
                        )


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


def insert_new_lineages(ncbi, df, tax_rank_filter, flowcell):
    """
    Insert any new lineages of species that are in our results into the database, and fetch all that we need for these
    results
    :param ncbi: The instance of the NCBITaxa class
    :param df: The data frame of our results
    :param tax_rank_filter: The order of our taxonomic ranks
    :param flowcell: The flowcell object
    :return lineages_Df: a dataframe of all the lineages
    """
    # create new defaultDict which creates a default dict for missing values
    lineage_dict = defaultdict(lambda: defaultdict())
    # Get the current taxIDs in the Dataframe
    taxid_list = df["tax_id"].values
    # Get the taxIDs of lineages that are in the database, so we can use existing lineages or create new ones
    already_lineaged_tax_ids = LineageValue.objects.values_list("tax_id", flat=True)

    # Check the taxIDs we have in the dataframe, to see if they have matching lineages already.
    # If not the set subtraction means that unique, lineage-less taxIDs end up in the not_prev_lineaged dict
    not_prev_lineaged = list(set(taxid_list) - set(already_lineaged_tax_ids))

    # Get the lineages for these taxIDs in the new centrifuge output,
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
        lineage_dict[key] = {rank: taxid_species_lookup_dict[tax_id]
                             for tax_id, rank in lineage_ranked.items()
                             if rank in tax_rank_filter or key == tax_id}

    # create dataframe, the .T means that the dataframe is transposed, so dict key is the series title
    lineages_df = pd.DataFrame.from_records(lineage_dict).T
    # rename the no rank column to strain
    lineages_df = lineages_df.rename(columns={"no rank": "strain"})
    # if these series don't exist add them and populate with numpy NaN so there's no key errors
    for series in ["subspecies", "strain", "subStrainSpecies"]:
        if series not in lineages_df.keys():
            lineages_df[series] = np.NaN

    # some strains are strains of a subspecies. But ncbi only returns the terminal (leaf) which is the
    # strain and not the subspecies above it so this function splits out the subspecies and adds it to a new
    # series
    new_lineages_for_saving_df = lineages_df[lineages_df.index.isin(not_prev_lineaged)]
    # if the df isn't empty
    if not new_lineages_for_saving_df.empty:
        # iterate over lineage dataframe to save to database #TODO rewrite into using mysql
        # set the index name so upon reset the added series has a name
        new_lineages_for_saving_df.index.name = "tax_id"

        new_lineages_for_saving_df.reset_index(inplace=True)
        # list to contain the bulk insert object
        lineages_insert_list = new_lineages_for_saving_df.apply(create_lineage_models,
                                                                axis=1)
        # Bulk create
        LineageValue.objects.bulk_create(list(lineages_insert_list.values))

    logger.info("Flowcell id: {} - Determining Subspecies".format(flowcell.id))
    # create new additional subspecies column, from strains column which has subspecies above it in taxa level
    lineages_df["subStrainSpecies"] = lineages_df["strain"].map(convert_species_to_subspecies)
    # merge new subspecies column with existing column
    lineages_df["subspecies"].fillna(lineages_df["subStrainSpecies"], inplace=True)

    # delete the new subspecies column
    delete_series(["subStrainSpecies"], lineages_df)

    return lineages_df


def calculate_sankey_values(lineages_df, df, flowcell, tax_rank_filter, job_master):
    """
    Calculate the values for the sankey diagram and save them into the database.
    :param lineages_df: A dataframe of lineages created from the current centrifuge output
    :param df: The main dataframe of centrifuge results
    :param flowcell: The flowcell for these reads
    :param tax_rank_filter: The taxonomic ranks in descending order
    :param job_master: The jobmaster object for this task
    :return:
    """
    # ##############  Storing snakey links and nodes ################

    logger.info("Flowcell id: {} - Starting Sankey Calculations".format(flowcell.id))
    # pd.merge outer join to combine the rows on tax_id, one row for each tax_id
    merged_lineages_df = pd.merge(lineages_df, df, how="outer", left_index=True, right_on="tax_id")
    # set index to tax_id
    merged_lineages_df.set_index("tax_id", inplace=True)

    # The columns we actually need
    columns = ["num_matches", "sum_unique", "barcode_name", "superkingdom", "phylum", "class", "order", "family",
               "genus", "species"]

    # Subset the dataframe into those columns
    sankey_lineages_df = merged_lineages_df[columns]

    # Backfill the dataframe nas with the lowest common entry, if family is missing, and genus
    #  isn't, family NaN becomes the genus entry
    sankey_lineages_df = sankey_lineages_df.fillna(axis=1, method="bfill")

    # initialise the series that will become our sankey pandas DataFrame
    source = pd.Series()
    target = pd.Series()
    value = pd.Series()
    tax_id = pd.Series()
    barcode = pd.Series()
    target_tax_level = pd.Series()
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

    # Create a DataFrame of links
    source_target_df = pd.concat([source, target, value, tax_id, target_tax_level, barcode], axis=1)

    # Rename the columns
    source_target_df.columns = ["source", "target", "value", "tax_id", "target_tax_level", "barcode_name"]

    # Create a series of Job_master objects through broadcasting
    source_target_df["job_master"] = job_master

    logger.info("Flowcell id: {} - Barcodes are:".format(flowcell.id))
    logger.info(source_target_df["barcode_name"].unique())

    # TODO DO sorting magic before the database deposition
    # Create a dataframe of the previous results for this metagenomics task
    previous_link_df = pd.DataFrame(
        list(SankeyLink.objects.filter(task__id=job_master.id).values())
    )

    # If there no results, i.e this is the first iteration
    if previous_link_df.empty:
        # The entire dataframe needs to be added into the database
        to_create_sankey_df = source_target_df
    # if there are results these need updating where we have the same link
    else:
        tax_id_list = list(previous_link_df["tax_id"])

        barcodes_list = list(previous_link_df["barcode_name"])
        # list of tuples
        tupley_tuple = list(zip(tax_id_list, barcodes_list))
        # Mask whether there is a link for this barcode and tax_id
        prev_links_mask = source_target_df[['tax_id', 'barcode_name']].agg(tuple, 1).isin(tupley_tuple)
        # subset where there is so these need dootdating
        to_create_sankey_df = source_target_df[~prev_links_mask]

    # ## ### Create the SankeyLinks object #### # One for each link #### #
    to_create_sankey_df.reset_index(inplace=True)

    # Remove any in this iteration where the value is 4 or less, as these will not be in the top 30
    to_create_sankey_df = to_create_sankey_df[to_create_sankey_df["value"] > 2]
    # Drop duplicates as we only need to create one
    to_create_sankey_df = to_create_sankey_df.drop_duplicates(subset=["tax_id", "target_tax_level"])

    # Create a second dataframe for values that need updating
    to_update_sankey_df = pd.DataFrame()

    # If there are values that need updating
    if not previous_link_df.empty:
        # Merge inner so only rows that are in both dataframe are kept, these are the rows that need updating
        # The target taxa level makes sure the rows intersect to the correct rank
        to_update_sankey_df = pd.merge(source_target_df, previous_link_df, how="inner", on=["tax_id", "barcode_name",
                                                                                            "target_tax_level"])

    # If there are values that need creating
    if not to_create_sankey_df.empty:
        # Perform the apply to create the series to create
        sankey_models_for_bulk_insert = to_create_sankey_df.apply(create_sankeylink_models, axis=1)

        # Bulk create the objects
        SankeyLink.objects.bulk_create(list(sankey_models_for_bulk_insert.values))

    # If there are values that need updating
    if not to_update_sankey_df.empty:
        # Sum the two values from the merged previous and current dataframes together
        to_update_sankey_df["updated_value"] = to_update_sankey_df["value_x"] + to_update_sankey_df["value_y"]
        logger.info("Flowcell id: {} - Updating Sankey Links - This will take a minute...".format(
            flowcell.id)
        )
        # TODO slowest point in code
        # Apply the update function to the dataframe
        to_update_sankey_df.apply(update_sankeylink_values, args=(
            job_master.id,), axis=1)
        logger.info("Flowcell id: {} - Finished updating sankey links".format(flowcell.id))


def insert_or_update_num_mapped(row, species_to_pk):
    """
    Insert the model object for the reads that mapped to their target reference, outside target regions
    :param row: Dataframe row
    :param species_to_pk: The primary key for each species MappingResults database row in a dict
    :return:
    """
    # TODO a lot of database calls, may need changing
    pk = species_to_pk[row["name"]]
    obj, created = MappingResultsBarcoded.objects.get_or_create(
        tax_id=row["tax_id"],
        mapping_result=pk,
        barcode_name=row["barcode_name"],
        species=row["name"]
    )
    obj.num_mapped += row["num_mapped"]
    obj.save()


def insert_or_update_num_mapped_barcoded(barcode_group_df, species_to_pk):
    """
    Calculate the num_mapped for reads that mapped to their reference outside target regions
    :param barcode_group_df: The dataframe for each group of barcodes
    :param species_to_pk: The primary key of the species entry in Cartography mapped, to use as Foreign Key
    :return:
    """
    name_map_gb = barcode_group_df.groupby("name")

    barcode_group_df.set_index("name", inplace=True)

    barcode_group_df["num_mapped"] = name_map_gb.size()

    barcode_group_df.reset_index(inplace=True)

    barcode_group_df.apply(insert_or_update_num_mapped, args=(species_to_pk,), axis=1)

    return


def plasmid_mapping(row, species, reference_location, fastq, flowcell):
    """
    Map the reads for groups that have danger regions on plasmids
    :param row: The row of the target dataframe for species with plasmid danger regions
    :param species: The species name
    :param reference_location: The location of the reference files on the system
    :param fastq: The reads sequences
    :param flowcell: The flowcell for logging by flowcell_id
    :return:
    """
    species = species.replace(" ", "_")

    reference_name = species + "_" + row["name"]

    references = ReferenceInfo.objects.get(name=reference_name)

    minimap2_ref_path = reference_location + references.filename

    map_cmd = '{} -x map-ont -t 4 --secondary=no {} -'.format("minimap2", minimap2_ref_path)

    # Execute minimap2
    out, err = subprocess.Popen(map_cmd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE) \
        .communicate(input=str.encode(fastq))
    # Get the output of minimap2
    map_output = out.decode()

    if not map_output:
        logger.info("Flowcell id: {} - No mappings for plasmid {} on species {}".format(flowcell.id, row["name"],
                                                                                        species))
        return
    else:
        plasmid_map_df = pd.read_csv(StringIO(map_output), sep="\t", header=None)

        columns = ["read_id", "query_seq_len", "query_start", "query_end", "rel_strand", "target_seq_name",
                   "target_seq_len", "target_start", "target_end", "num_matching_bases", "num_matching_bases_gaps",
                   "mapping_qual", "type_align", "number_minimiser", "chaining_score", "chaining_score_2nd_chain",
                   "random"]
        # set the column headers to above
        plasmid_map_df.columns = columns
        # Drop unnecessary columns
        plasmid_map_df.drop(columns=["number_minimiser", "chaining_score", "chaining_score_2nd_chain", "random"],
                            inplace=True)
        # TODO NOTE that atm, we're just counting any map that reads to a plasmid as a red read

        logger.info(
            "Flowcell id: {} - This many reads mapped evilly on this reference {} for species {}"
            .format(flowcell.id, plasmid_map_df.shape[0], species)
        )
        # return as tuple to iterate over
        return plasmid_map_df


def map_all_the_groups(target_species_group_df, group_name, reference_location, flowcell, gff3_df, targets_results_df,
                       task, species_to_pk):
    """
    Map the reads from the target data frames, after they've been grouped by species
    :param target_species_group_df: A data frame that contains reads from only one species
    :param group_name: The name of the species in this group
    :param reference_location: Where we have stored our reference files for minimap2
    :param flowcell: The flowcell model from where these reads came from
    :param gff3_df: The data frame containing the target regions
    :param targets_results_df: The data frame containing the new targets mapping against reference
    :param task: The task object
    :param species_to_pk: The primary key to species dict lookup
    :return red_df: Any reads that map to plasmids
    """
    # The species as identified by centrifuge for this group of reads

    logger.info("Flowcell id: {} - species is {}".format(flowcell.id, group_name))
    # The reference file for this species
    species_regions_df = gff3_df[gff3_df["species"] == group_name]

    species = group_name.replace(" ", "_")

    try:
        references = ReferenceInfo.objects.get(name=species)

    except ObjectDoesNotExist:

        logger.info("Flowcell id: {} - No reference found for species {}!".format(flowcell.id, species))

        return

    # the minimap2 reference fasta
    minimap2_reference_path = reference_location + references.filename
    # The read sequences for the reads we want to map
    reads = FastqRead.objects.filter(run__flowcell_id=flowcell.id, read_id__in=target_species_group_df["read_id"])
    # The fasta sequence
    fastqs_list = reads.values_list('read_id', 'fastqreadextra__sequence')
    # Assemble the fastq into a string
    fastq = "".join([str(">" + c[0] + "\n" + c[1] + "\n")
                     for c in fastqs_list if c[1] is not None])
    # The command to execute minimap
    map_cmd = '{} -x map-ont -t 4 --secondary=no {} -'.format("minimap2", minimap2_reference_path)
    # Execute minimap2
    out, err = subprocess.Popen(map_cmd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE) \
        .communicate(input=str.encode(fastq))
    # Get the output of minimap2
    map_output = out.decode()

    logger.info("Flowcell id: {} - map out {} ".format(flowcell.id, map_output))

    plasmid_red_df = pd.DataFrame()

    plasmid_df = species_regions_df[species_regions_df["gff_line_type"] == "plasmid"]
    # If there is output from minimap2 create a dataframe
    if not plasmid_df.empty:
        logger.info("Flowcell id: {} - Mapping reads to plasmids for species {} ".format(flowcell.id, species))
        # TODO is this doubling the first line?
        plasmid_red_df = plasmid_df.apply(plasmid_mapping, args=(species, reference_location, fastq, flowcell),
                                          axis=1)

    if map_output:
        map_df = pd.read_csv(StringIO(map_output), sep="\t", header=None)
    # Else stop and return. There is no output if there are no mappings
    else:
        logger.info("Flowcell id: {} - No minimap output (No reads mapped) for species {}".format(flowcell.id, species))
        return
    # The columns headers
    columns = ["read_id", "query_seq_len", "query_start", "query_end", "rel_strand", "target_seq_name",
               "target_seq_len", "target_start", "target_end", "num_matching_bases", "num_matching_bases_gaps",
               "mapping_qual", "type_align", "number_minimiser", "chaining_score", "chaining_score_2nd_chain",
               "random"]
    # set the column headers to above
    map_df.columns = columns
    # Drop unecessary columns
    map_df.drop(columns=["number_minimiser", "chaining_score", "chaining_score_2nd_chain", "random"],
                inplace=True)

    boolean_df = species_regions_df.apply(falls_in_region, args=(map_df,), axis=1)
    # Get whether it's inside the boundaries or not
    boolean_df = boolean_df.any()
    # Red reads, reads that fall within boundaries
    red_df = map_df[boolean_df]
    if not plasmid_red_df.empty:
        red_df = red_df.append(plasmid_red_df)

    logger.info(
        "Flowcell id: {} - This many reads mapped evilly on this reference {} for species {}".format(
            flowcell.id,
            red_df.shape[0],
            species)
    )
    # Non red reads, reads that mapped elsewhere on the reference
    map_df = map_df[~boolean_df]

    logger.info(
        "Flowcell id: {} - This many reads mapped elsewhere on this reference {}".format(flowcell.id,
                                                                                         map_df.shape[0]))
    # get the info on the mapped reads from the original dataframe
    # TODO this is where we sort the barcodes
    update_num_mapped_df = pd.merge(targets_results_df, map_df, how="inner", left_on="read_id", right_on="read_id")
    # Group by name, to get size of group, so number of reads matching
    gb_species_mapped_non_read = update_num_mapped_df.groupby("name")
    # Set the index as name
    update_num_mapped_df.set_index("name", inplace=True)
    # Get the number of reads that mapped to each species
    update_num_mapped_df["num_mapped"] = gb_species_mapped_non_read.size()
    # Drop duplicates so unique in each species
    update_num_mapped_df = update_num_mapped_df[~update_num_mapped_df.index.duplicated(keep="first")]

    update_num_mapped_df.reset_index(inplace=True)
    # Apply to each row to make an update
    update_num_mapped_df.apply(update_mapped_non_dangerous, args=(task,), axis=1)

    # Update the barcoded mapped but none dangerous reads
    gb_barcode_mapped_non_red = update_num_mapped_df.groupby("barcode_name")

    for name, group in gb_barcode_mapped_non_red:
        insert_or_update_num_mapped_barcoded(group, species_to_pk)
    # gb_unred_mapped_barcode.apply(insert_or_update_num_mapped_barcoded, species_to_pk)

    logger.info("Flowcell id: {} - Finished updating mapped reads".format(flowcell.id))

    return red_df


def falls_in_region(row, map_df):
    """
    Does this reads mapping fall within any of the regions we have defined?
    :param row: The gff3 dataframe row
    :param map_df: The results of the mapping step in a pandas dataframe
    :return:
    """
    # TODO which region?
    # If the start is higher than the region start
    start_low = map_df["target_start"] > row["start"]
    # If the start is lower than the region end
    start_high = map_df["target_start"] < row["end"]
    # If the end is higher than the region end
    end_low = map_df["target_end"] > row["start"]
    # If the end is lower than the region end
    end_high = map_df["target_end"] < row["end"]
    # Concat the series together, create "truth" dataframe
    bool_df = pd.concat([start_low, start_high, end_low, end_high], axis=1, ignore_index=True)
    # Keep where the start is in the region,
    bool_df["keep"] = np.where((bool_df[0] & bool_df[1]) | (bool_df[2] & bool_df[3]), True, False)
    return bool_df["keep"]


def update_mapped_red_values(row, bulk_insert_red_id, mapped_targets_dict, task):
    """
    Update the values that we have in the database for each species after getting results for this iteration
    :param row: The results_df row
    :param bulk_insert_red_id: The list to store the red reads database object in
    :param mapped_targets_dict: Whether we have update this entry or not as each species row only differs by read_id
    :param task: The task model object
    :return:
    """
    if row["name"] not in mapped_targets_dict:
        mapped = MappingResult.objects.get(species=row["name"],
                                           task=task)
        mapped.red_reads = row["summed_red_reads"]
        mapped.red_sum_unique = row["summed_sum_unique"]

        mapped.save()

        mapped_targets_dict[row["name"]] = mapped.id

        target_read_id = TargetMappedReadId(read_id=row["read_id"], mapping_result=mapped.id)

        bulk_insert_red_id.append(target_read_id)
    else:
        mapped_result = TargetMappedReadId(read_id=row["read_id"], mapping_result_id=mapped_targets_dict[row["name"]])

        bulk_insert_red_id.append(mapped_result)

    return bulk_insert_red_id


def update_mapped_non_dangerous(row, task):
    """
    Update number of mapped reads that don't map to danger but do map to reference
    :param row: DataFrame row
    :param task: The related task for this analysis
    :return:
    """
    mapped_result = MappingResult.objects.get(task=task, species=row["name"])
    mapped_result.num_mapped += row["num_mapped"]
    mapped_result.save()


def create_barcode_mapped_dangerous(row, species_to_pk):
    """
    Insert the updated values into the database
    :param row: The dataframe row
    :param species_to_pk: A dict with lookup for species to their database row primary key, used to Foreign key
    to Cartography mapped
    :return:
    """
    pk = species_to_pk[row["name"]]

    obj, created = MappingResultsBarcoded.objects.get_or_create(
        tax_id=row["tax_id"],
        mapping_result=pk,
        barcode_name=row["barcode_name"],
        species=row["name"]
    )

    obj.red_reads += row["barcode_red_reads"]

    obj.red_sum_unique += row["barcode_summed_red"]
    obj.save()


def update_mapped_dangerous_barcode(mapped_dangerous_barcode_group, species_to_pk):
    """
    Update or create the red reads values for Each Barcode
    :param mapped_dangerous_barcode_group: Groups of dangerous reads mapping to t
    :param species_to_pk: A dict with lookup for species to their database row primary key, used to Foreign key
    to Cartography mapped
    :return:
    """
    name_gb = mapped_dangerous_barcode_group.groupby("name")
    mapped_dangerous_barcode_group.set_index("name", inplace=True)
    mapped_dangerous_barcode_group["barcode_red_reads"] = name_gb.size()
    mapped_dangerous_barcode_group["barcode_summed_red"] = name_gb["unique"].sum()
    mapped_dangerous_barcode_group.reset_index(inplace=True)
    mapped_dangerous_barcode_group.drop_duplicates(subset=["barcode_name", "name"], inplace=True)
    mapped_dangerous_barcode_group.apply(create_barcode_mapped_dangerous, args=(species_to_pk,), axis=1)


def map_the_reads(name_df, task, flowcell):
    """
     Run the mapping step
    :return: The
    """
    # Targets_df
    logger.info("Flowcell id: {} - Mapping the target reads".format(flowcell.id))

    targets_df = name_df
    # TODO currently hardcoded below
    target_set = "starting_defaults"
    # if there are no targets identified by centrifuge in this iteration
    if targets_df.empty:

        logger.info("Flowcell id: {} - No targets in this batch of reads".format(flowcell.id,
                                                                                 targets_df.shape))
        # Get the targets set in the gff files and uploaded
        target_regions_df = pd.DataFrame(list(MappingTarget.objects.filter(target_set=target_set).values()))
        # If there are no target regions found for the set name
        if target_regions_df.empty:
            logger.info("Flowcell id: {} - No set of target regions found by that set name"
                        .format(flowcell.id, ))
            return
        # Get np array of the target species
        # TODO split into function so code is DRYER
        target_species = target_regions_df["species"].unique()
        # Get np array of their tax_ids
        target_tax_id = target_regions_df["tax_id"].unique()
        # Combine into tuple, one for each
        target_tuple = list(zip(target_species, target_tax_id))
        # create one MappingResults entry for each target species
        for target in target_tuple:
            obj, created = MappingResult.objects.get_or_create(task=task,
                                                               species=target[0],
                                                               tax_id=target[1],
                                                               defaults={"red_reads": 0,
                                                                         "num_mapped": 0,
                                                                         "red_sum_unique": 0}
                                                               )
            if created:
                logger.info("Flowcell id: {} - CM object created for {}".format(flowcell.id, obj.species))
        return

    # Reset the index
    targets_df.reset_index(inplace=True)

    logger.info("Flowcell id: {} - targets_df shape is {}".format(flowcell.id, targets_df.shape))
    # Location of the reference file
    reference_location = getattr(settings, "REFERENCE_LOCATION", None)

    queryset = MappingResult.objects.filter(task=task)
    # Create a lookup of pks to species
    species_to_pk = {query.species: query for query in queryset}
    # Get one group for each name
    gb = targets_df.groupby(["name"], as_index=False)
    # Get the targets set in the gff files and uploaded
    target_regions_df = pd.DataFrame(list(MappingTarget.objects.filter(target_set=target_set).values()))
    # Get np array of the target species
    # TODO split into function so code is DRYER
    target_species = target_regions_df["species"].unique()
    # Get np array of their tax_ids
    target_tax_id = target_regions_df["tax_id"].unique()
    # Combine into tuple, one for each
    target_tuple = list(zip(target_species, target_tax_id))

    # create one MappingResults entry for each target species
    for target in target_tuple:
        obj, created = MappingResult.objects.get_or_create(
            task=task,
            species=target[0],
            tax_id=target[1],
            defaults={"red_reads": 0,
                      "num_mapped": 0,
                      "red_sum_unique": 0}
        )
        if created:
            logger.info("Flowcell id: {} - CM object created for {}".format(flowcell.id, obj.species))
    # Initialise red_df for red alert target reads
    red_reads_df = pd.DataFrame()
    for name, group in gb:
        red_reads_df = red_reads_df.append(map_all_the_groups(group, name, reference_location, flowcell,
                                                              target_regions_df, targets_df, task,
                                                              species_to_pk)
                                           )

    # If none of the reads mapped to inside the groups
    if red_reads_df.empty:
        logger.info("Flowcell id: {} - No reads mapped to dangerous areas! Better luck next time."
                    .format(flowcell.id))
        return

    red_reads_df.set_index(["read_id"], inplace=True)
    # results df contains all details on red reads
    results_df = pd.merge(targets_df, red_reads_df, how="inner", left_on="read_id", right_index=True)
    # If there are no reads that are red reads
    mapping_gb = results_df.groupby(["name"])

    results_df.set_index(["name"], inplace=True)
    # Get how many reads are mapped to each genome dangerously
    results_df["red_num_matches"] = mapping_gb.size()
    # Get the sum of the unique reads
    results_df["red_sum_unique"] = mapping_gb["unique"].sum()

    results_df.reset_index(inplace=True)
    logger.info("Flowcell id: {} - The results of the mapping stage are {}".format(flowcell.id, results_df))
    logger.info("Flowcell id: {} - The results keys are {}".format(flowcell.id, results_df.keys()))

    # List to store the read read ids
    bulk_insert_red_read_ids_from_update = []

    mapped_target_dict = {}
    # Get the previous values for the mappings
    prev_df = pd.DataFrame(list(MappingResult.objects.filter(task=task).values()))
    prev_df.set_index(["tax_id"], inplace=True)
    logger.info("Flowcell id: {} - The previous results are {}".format(flowcell.id, prev_df.head()))
    logger.info("Flowcell id: {} - The previous results keys are {}".format(flowcell.id, prev_df.keys()))

    if not results_df.empty:
        results_df.set_index(["tax_id"], inplace=True)
        # Add the previous number of red reads to the results df
        results_df["to_add_red_reads"] = prev_df["red_reads"]
        # Add the previous sum unique to the results dataframe
        results_df["to_add_red_sum_unique"] = prev_df["red_sum_unique"]
        # Combine the previous and new number matches
        results_df["summed_red_reads"] = results_df["red_num_matches"] + results_df["to_add_red_reads"]
        # Combine the previous sum of the unique reads
        results_df["summed_sum_unique"] = results_df["red_sum_unique"] + results_df["to_add_red_sum_unique"]
        logger.info("Flowcell id: {} -  {}".format(flowcell.id, prev_df.keys()))
        logger.info("Flowcell id: {} -  {}".format(flowcell.id, prev_df.keys()))
        results_df.reset_index(inplace=True)
        # Apply the update mapped values function
        results_df.apply(update_mapped_red_values, args=(bulk_insert_red_read_ids_from_update,
                                                         mapped_target_dict,
                                                         task), axis=1)
        # BUlk create the red read ids
        try:
            TargetMappedReadId.objects.bulk_create(bulk_insert_red_read_ids_from_update)
        except IntegrityError:
            logger.info("Flowcell id: {} - duplicate red reads!".format(flowcell.id))
        # Apply the calculation and update to the dataframe
        barcode_mapped_gb = results_df.groupby("barcode_name")
        # Apply to each barcode group
        barcode_mapped_gb.apply(update_mapped_dangerous_barcode, species_to_pk)


def run_centrifuge(flowcell_job_id):
    """

    Returns nothing.
    -------
     Run the centrifuge command, and open the subprocess to run it. Keep the process open and pipe the reads into
     it Write the reads incrementally to the processes stdin, capture the stdout and parse it into a useful format.
     Once parsed it's deposited into the database, updating any existing records. The organising key is a uuid,
     stored in self.foreign key.
    """

    job_master = JobMaster.objects.get(pk=flowcell_job_id)
    job_master.running = True
    job_master.save()

    flowcell = job_master.flowcell

    logger.info('Flowcell id: {} - Running centrifuge on flowcell {}'.format(flowcell.id, flowcell.name))
    logger.info('Flowcell id: {} - job_master_id {}'.format(flowcell.id, job_master.id))
    logger.info('Flowcell id: {} - last_read {}'.format(flowcell.id, job_master.last_read))
    logger.info('Flowcell id: {} - read_count {}'.format(flowcell.id, job_master.read_count))

    centrifuge_path = get_env_variable("MT_CENTRIFUGE")
    index_path = get_env_variable("MT_CENTRIFUGE_INDEX")
    cmd = "perl " + centrifuge_path + " -f --mm -x " + index_path + " -"

    # Counter for total centOut
    total_centrifuge_output = 0

    # Instance of the Ncbi taxa class, for taxonomic id manipulation
    ncbi = NCBITaxa()

    # Get the task record from the JobMaster table in the database. Used to separate the results inserted,
    # and tell the Javascript and RunMonitor the task is running and when it is complete

    fastqs = FastqRead.objects.filter(run__flowcell=flowcell, id__gt=int(job_master.last_read)
                                      ).order_by('id')[:2000]

    logger.info('Flowcell id: {} - number of reads found {}'.format(flowcell.id, fastqs.count()))

    try:
        metadata = Metadata.objects.get(task=job_master)

    except (Metadata.DoesNotExist, Metadata.MultipleObjectsReturned):

        metadata = Metadata(

            task=job_master
        )

        metadata.save()

    fastqs_list = fastqs.values_list('read_id', 'barcode_name', 'fastqreadextra__sequence')

    fastqs_data = "".join([str(">read_id=" + c[0] + ",barcode=" + c[1] + "\n" + c[2] + "\n") for c in fastqs_list
                           if c[2] is not None])

    # Write the generated fastq file to stdin, passing it to the command
    logger.info("Flowcell id: {} - Loading index and Centrifuging".format(flowcell.id))

    # Use Popen to run the centrifuge command
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE) \
        .communicate(input=str.encode(fastqs_data))

    # the number of centrifuge lines
    centrifuge_output = out.decode()

    # total number of lines of centrifuge output dealt with
    total_centrifuge_output += centrifuge_output.count("\n") - 1

    logger.info("Flowcell id: {} - number of centrifuge output lines is {}".format(flowcell.id, total_centrifuge_output))

    # output fields is the column headers for the pandas data frame
    output_fields = ["readID", "seqID", "taxID", "numMatches"]

    # create the DataFrame from the output
    df = pd.read_csv(StringIO(centrifuge_output), sep="\t", usecols=output_fields)

    if df.empty:
        logger.info("Flowcell id: {} - No reads found or no centrifuge output".format(flowcell.id))
        return

    # create DataFrame column unique, where if the number of matches is 1, unique is 1, if matches is more than
    # 1 unique is 0, using numpy condition
    df["unique"] = np.where(
        df["numMatches"] == 1,  # condition
        1,  # True
        0  # False
    )

    df = pd.concat([df, df["readID"].str.split(",").str[1].str.split("=").str[1]], axis=1)

    df.columns = ["readID", "seqID", "tax_id", "numMatches", "unique", "barcode_name"]

    df = pd.concat([df, df["readID"].str.split(",").str[0].str.split("=").str[1]], axis=1)

    df.columns = ["readID", "seqID", "tax_id", "numMatches", "unique", "barcode_name", "read_id"]

    # get np array of the taxids
    taxid_list = np.unique(df["tax_id"].values)

    # use the ncbi thing to get the species names for the tax ids in the tax ids list
    taxid2name_df = pd.DataFrame.from_dict(ncbi.get_taxid_translator(taxid_list), orient="index",
                                           columns=["name"])

    # insert the taxid2name dict items into the dataframe name columns
    df = pd.merge(df, taxid2name_df, how="outer", left_on="tax_id", right_index=True)

    # any NaNs replace with Unclassified
    df["name"].fillna("Unclassified", inplace=True)

    #  set the index to be taxID
    df = df.set_index("tax_id")
    # The barcoding
    gb_bc = df.groupby(["barcode_name"])

    # create lots of mini data frames grouping by taxID
    gb = df.groupby(level="tax_id")

    # create a sumUnique column in the DataFrame, where the value
    # is the sum of the unique column in the corresponding
    # grouped_by table
    df["sum_unique"] = gb["unique"].sum()

    # create a number of reads column in the data frame by getting
    # the number of rows in the corresponding grouped_by
    # table
    df["num_matches"] = gb.size()

    # TODO future mapping target

    # Unique targets in this set of regions
    temp_targets = MappingTarget.objects.filter(target_set="starting_defaults") \
        .values_list("species", flat=True).distinct()

    logger.info("Flowcell id: {} - Targets are {}".format(flowcell.id, temp_targets))

    name_df = df[df["name"].isin(temp_targets)]

    logger.info("Flowcell id: {} - The dataframe target dataframe is {}".format(flowcell.id, name_df))

    map_the_reads(name_df, job_master, flowcell)

    barcode_df = pd.DataFrame()

    # TODO vectorise these bad boys
    logger.info("Flowcell id: {} - The dataframe shape is {}".format(flowcell.id, df.shape))

    unique_barcode = set([c[1] for c in fastqs_list])

    logger.info("Flowcell id: {} - These are the barcodes in this set {}".format(flowcell.id, unique_barcode))

    barcode_df = gb_bc.apply(calculate_barcoded_values, barcode_df)

    barcode_df.reset_index(inplace=True)

    barcode_df.rename(columns={"unique": "sum_unique", 0: "num_matches"}, inplace=True)

    barcode_df.set_index("tax_id", inplace=True)

    # delete these columns, no longer needed
    df.drop(columns=["readID", "seqID", "numMatches", "unique", "barcode_name", "read_id"], inplace=True)

    # remove duplicate rows in the data frame,so we only have one entry for each species
    df.drop_duplicates(keep="first", inplace=True)

    # filter out no ranks, like root and cellular organism
    tax_rank_filter = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

    # drop all duplicated lines to keep only one for each entry, Atomicity
    df.reset_index(inplace=True)

    df = df.drop_duplicates(keep="first")
    # Make this the results for all reads in the centrifuge output
    df["barcode_name"] = "All reads"
    # query the index get all the objects in the database for this analysis at this point
    queryset = CentrifugeOutput.objects.filter(task=job_master)
    # Get the previous tax ids in this analysis
    prev_df_tax_ids = list(queryset.values_list("tax_id", flat=True))

    # ###### BULK CREATE CENTOUTPUT OBJECTS ########
    # Give each row a flowcell object cell, as part of the flowcell_id series
    df["flowcell"] = flowcell

    # same with task from Jobmaster
    df["task"] = job_master

    # Split out the rows that have a tax_id that doesn't already have results
    cent_to_create_df = df[~df["tax_id"].isin(prev_df_tax_ids)]

    # apply row wise to append model representation objects into a series
    centrifuge_create_series = cent_to_create_df.apply(create_centrifuge_models, axis=1)

    logger.info("Flowcell id: {} - Bulk creating CentOutput objects".format(flowcell.id))
    # Bulk create the objects
    CentrifugeOutput.objects.bulk_create(list(centrifuge_create_series.values))

    df.set_index("tax_id", inplace=True)
    # Add the barcode dataframe onto the dataframe, so now we have all reads and barcodes
    df = df.append(barcode_df, sort=True)

    df.reset_index(inplace=True)
    # Get all the results for the barcoded entries, so each species will have multiple entries under different barcodes
    # TODO Split back into one table in database, barcoded results
    bar_queryset = CentrifugeOutputBarcoded.objects.filter(output__task__id=job_master.id).values()
    # If there are results, we need to update existing results database rows and create new ones for non extant rows
    if bar_queryset:
        logger.info("Flowcell id: {} - Previous CentOutput data found".format(flowcell.id))
        # Create a pandas dataframe
        to_create_bar_df = pd.DataFrame(list(bar_queryset))
        # These results will be removed from the dataframe
        to_create_bar_df["temp"] = "Y"
        # The original dataframe results will be kept
        df["temp"] = "N"
        # Append the two dataframes together
        to_create_bar_df = df.append(to_create_bar_df, sort=True)
        # separate out lines with duplicated tax_id and barcodes, these are lines with values that need updating in the
        #  database, keep both lines
        to_update_bar_df = to_create_bar_df[to_create_bar_df.duplicated(subset=["tax_id", "barcode_name"],
                                                                        keep=False)]
        # set multi index to sort by tax_ids into each barcode
        to_update_bar_df.set_index(["barcode_name", "tax_id"], inplace=True)
        # Sum the sum unique for each species in each barcode
        to_update_bar_df["updated_sum_unique"] = to_update_bar_df.groupby(["barcode_name",
                                                                           "tax_id"])["sum_unique"].sum()
        # Sum the number of matches as above. Same as fr all reads calculations
        to_update_bar_df["updated_num_matches"] = to_update_bar_df.groupby(["barcode_name",
                                                                            "tax_id"])["num_matches"].sum()
        # Remove duplicated lines
        to_update_bar_df = to_update_bar_df[to_update_bar_df.index.duplicated(keep="first")]

        to_update_bar_df.reset_index(inplace=True)
        # ###### Update existing barcode output ######
        # Apply to update the existing entries in the databases
        to_update_bar_df.apply(update_centrifuge_barcoded_values, args=(job_master.id,),
                               axis=1)
        # Drop duplicates of barcodes under each tax_id so we only create one entry to update
        to_create_bar_df.drop_duplicates(subset=["tax_id", "barcode_name"], inplace=True, keep=False)
        # Remove original df entries, so just barcoded data, not "All reads"
        to_create_bar_df = to_create_bar_df[to_create_bar_df["temp"] == "N"]
        # Create the model objects, returns a series
        centrifuge_models_series = to_create_bar_df.apply(create_centrifuge_barcoded_models, args=(job_master,),
                                                      axis=1)
        # Bulk create the value
        CentrifugeOutputBarcoded.objects.bulk_create(list(centrifuge_models_series.values))
        logger.info("Flowcell id: {} - Bulk creating CentOutputBarcoded objects".format(
            flowcell.id)
        )
        # remove unnecessary temp series
        df.drop(columns=["temp"], inplace=True)

    else:
        # apply to create the models into a series
        centrifuge_models_series = df.apply(create_centrifuge_barcoded_models, args=(job_master,),
                                        axis=1)
        # Bulk create
        CentrifugeOutputBarcoded.objects.bulk_create(list(centrifuge_models_series.values))

        logger.info("Flowcell id: {} - Bulk creating CentOutputBarcoded objects".format(
            flowcell.id)
        )

    lineages_df = insert_new_lineages(ncbi, df, tax_rank_filter, flowcell)
    calculate_sankey_values(lineages_df, df, flowcell, tax_rank_filter, job_master)

    # ##################### Back to the non sankey calculations ####################
    # Update the metadata for the headers
    logger.info("Flowcell id: {} - Finished all reads in database".format(flowcell.id))

    # Get the start time #TODO remove the unused fields
    start_time = metadata.start_time.replace(tzinfo=None)
    end_time = timezone.now().replace(tzinfo=None)

    metadata.finish_time = str(end_time - start_time)
    metadata.save()

    # Update the jobmaster object fields that are relevant
    job_master = JobMaster.objects.get(pk=job_master.id)
    job_master.running = False

    if fastqs.count() > 0:
        job_master.last_read = fastqs[499].id

    job_master.read_count = job_master.read_count + fastqs.count()
    job_master.save()

    logger.info("Flowcell id: {} - New last_read_id is - {}".format(flowcell.id, job_master.last_read))
    logger.info("Flowcell id: {} - Total CentOut lines - {}".format(flowcell.id, total_centrifuge_output))
    logger.info("Flowcell id: {} - Finished!".format(flowcell.id))
