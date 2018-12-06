"""Centrifuge.py
"""
import os
import subprocess
from collections import defaultdict
from io import StringIO
import numpy as np
import pandas as pd
from celery.utils.log import get_task_logger
from django.db.models import ObjectDoesNotExist
from django.utils import timezone
from ete3 import NCBITaxa
from centrifuge.models import CentrifugeOutput, LineageValue, Metadata, SankeyLink, \
    MappingResult, MappingTarget, DonutData  # , TargetMappedReadId
from jobs.models import JobMaster
from minotourapp.utils import get_env_variable
from reads.models import FastqRead
from reference.models import ReferenceInfo
from django.conf import settings
# from django.db.utils import IntegrityError
from django.db.models import Sum

pd.options.mode.chained_assignment = None
logger = get_task_logger(__name__)


# TODO del unused df and things

def calculate_barcoded_values(barcode_group_df, barcode_df, classified_per_barcode):
    """
    Return the main cent output data frame with the barcode only result concatenated on
    :param barcode_group_df: barcode data frame group data frame, one for each barcode
    :param barcode_df: the main results data frame that will store the ouput of this function
    :param classified_per_barcode: The number of classified reads in a barcode
    :return: A new df with the result for barcode appended
    """
    barcode = barcode_group_df["barcode_name"].unique()[0]

    reads_classed_in_barcode = classified_per_barcode[barcode]

    barcode_group_df_gb = barcode_group_df.groupby("tax_id")
    # calculate total unique matches for this species in this barcode as a series
    sum_unique = barcode_group_df_gb["unique"].sum()
    # calculate the total matches fot this species in this barcode as a series
    num_matches = barcode_group_df_gb.size()

    barcode_group_df.drop(columns=["readID", "seqID", "numMatches", "unique"], inplace=True)

    barcode_group_df.drop_duplicates(subset=["barcode_name", "name"], keep="first", inplace=True)
    # Calculate the proportion
    barcode_group_df["proportion_of_classified"] = barcode_group_df["num_matches"] \
        .div(reads_classed_in_barcode).mul(100).round(decimals=3)

    # combine the dataframes
    values_df = pd.concat([sum_unique, num_matches,
                           barcode_group_df["name"],
                           barcode_group_df["proportion_of_classified"]], axis=1)

    barcode_df = barcode_df.append(values_df, sort=True)

    return barcode_df


def create_centrifuge_models(row, classified_per_barcode):
    """
    Append a CentOutput object to a list for each row in the centrifuge output
    :param row: the row from the data frame
    :param classified_per_barcode: The number of reads classified for each barcode as a dictionary
    :return: The list of newly created objects

    """
    if row["proportion_of_classified"] == "Unclassified":
        row["proportion_of_classified"] = round(row["num_matches"] / classified_per_barcode[row["barcode_name"]], 3)
    return CentrifugeOutput(name=row["name"],
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
                            species=row["species"])


def update_centrifuge_output_values(row, flowcell_job_id):
    """
    Update existing CentrifugeOutputBarcoded objects in the database
    :param flowcell_job_id: The id for the job_master for this task run
    :param row: The data frame row
    :return: The list of newly created objects
    """
    CentrifugeOutput.objects.filter(tax_id=row["tax_id"],
                                    barcode_name=row["barcode_name"],
                                    task__id=flowcell_job_id
                                    ).update(
        num_matches=row["num_matches"],
        sum_unique=row["sum_unique"],
        proportion_of_classified=row["proportion_of_classified"]
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
                      barcode_name=row["barcode_name"],
                      path=row["path"])


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


def update_mapped_non_dangerous(row, task):
    """
    Update number of mapped reads that don't map to danger but do map to reference
    :param row: DataFrame row
    :param task: The related task for this analysis
    :return:
    """
    mapped_result = MappingResult.objects.get(task=task, species=row["name"], barcode_name=row["barcode_name"])
    mapped_result.num_mapped += row["num_mapped"]
    nm = mapped_result.num_mapped
    mapped_result.mapped_proportion_of_classified = round((nm / mapped_result.num_matches) * 100, 2)
    mapped_result.save()


def plasmid_mapping(row, species, fastq_list, flowcell):
    """
    Map the reads for groups that have danger regions on plasmids
    :param row: The row of the target dataframe for species with plasmid danger regions
    :param species: The species name
    :param fastq_list: The reads sequences
    :param flowcell: The flowcell for logging by flowcell_id
    :return plasmid_map_df: A list of dicts containing the information about the plasmid mappings
    """

    #
    # TODO Delete next lines
    #
    fastq_input = ''

    for fastq in fastq_list:

        if fastq[1]:

            temp = '>{}\n{}\n'.format(fastq[0], fastq[1])
            fastq_input += temp

    fastq = fastq_input

    #
    # End of deletion
    #

    reference_location = getattr(settings, "REFERENCE_LOCATION", None)
    minimap2 = getattr(settings, "MINIMAP2", None)

    species = species.replace(" ", "_")

    reference_name = species + "_" + row["name"]

    references = ReferenceInfo.objects.get(name=reference_name)

    minimap2_ref_path = os.path.join(reference_location, references.filename)

    map_cmd = '{} -x map-ont -t 4 --secondary=no {} -'.format(minimap2, minimap2_ref_path)
    logger.info(map_cmd)

    # Execute minimap2
    out, err = subprocess.Popen(map_cmd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                stderr=subprocess.PIPE).communicate(input=str.encode(fastq))
    if err:
        logger.info("Flowcell id: {} - Minimap error!! {}".format(flowcell.id, err))
    # Get the output of minimap2
    map_output = out.decode()
    logger.info("Flowcell id: {} - Plasmid mapping map output is {}".format(flowcell.id, map_output))
    columns = ["read_id", "query_seq_len", "query_start", "query_end", "rel_strand", "target_seq_name",
               "target_seq_len", "target_start", "target_end", "num_matching_bases", "num_matching_bases_gaps",
               "mapping_qual", "type_align", "number_minimiser", "chaining_score", "chaining_score_2nd_chain",
               "random"]
    if not map_output:
        logger.info("Flowcell id: {} - No mappings for plasmid {} on species {}".format(flowcell.id, row["name"],
                                                                                        species))
        return np.NaN
    else:
        plasmid_map_df = pd.read_csv(StringIO(map_output), sep="\t", header=None)

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
        logger.info(plasmid_map_df.head())

        return plasmid_map_df


def update_mapped_red_values(row, task, flowcell):
    """
    Update the values that we have in the database for each species after getting results for this iteration
    :param row: The results_df row
    :param task: The task model object
    :param flowcell: The flowcell model object
    :return:
    """
    logger.info("Flowcell id: {} - updating the number of mapped red reads for species {}"
                .format(flowcell.id, row["species"]))
    mapped = MappingResult.objects.get(species=row["species"],
                                       task=task, barcode_name=row["barcode_name"])
    nm = mapped.num_mapped
    mapped.red_reads += row["red_num_matches"]
    mapped.red_reads_proportion_of_classified = round((mapped.red_reads / nm) * 100, 2)
    mapped.save()


def map_all_the_groups(target_species_group_df, group_name,  flowcell, gff3_df, targets_results_df,
                       task, num_matches_target_barcoded_df):
    """
    Map the reads from the target data frames, after they've been grouped by species
    :param target_species_group_df: A data frame that contains reads from only one species
    :param group_name: The name of the species in this group
    :param flowcell: The flowcell model from where these reads came from
    :param gff3_df: The data frame containing the target regions
    :param targets_results_df: The data frame containing the new targets mapping against reference
    :param task: The task object
    :param num_matches_target_barcoded_df: The number of matches per barcode for each target species in a dataframe
    :return red_df: Any reads that map to plasmids and their information as a dataframe
    """
    # The species as identified by centrifuge for this group of reads

    logger.info("Flowcell id: {} - species is {}".format(flowcell.id, group_name))
    logger.info("Flowcell id: {} - The num_matches_df is ".format(flowcell.id, num_matches_target_barcoded_df))
    # The reference file for this species
    species_regions_df = gff3_df[gff3_df["species"] == group_name]

    species = group_name.replace(" ", "_")

    try:
        references = ReferenceInfo.objects.get(name=species)

    except ObjectDoesNotExist:

        logger.info("Flowcell id: {} - No reference found for species {}!".format(flowcell.id, species))

        return

    reference_location = getattr(settings, "REFERENCE_LOCATION", None)
    minimap2 = getattr(settings, "MINIMAP2", None)

    # the minimap2 reference fasta
    # minimap2_reference_path = os.path.join(reference_location, references.filename)
    minimap2_reference_path = os.path.join(reference_location, references.filename)
    # The read sequences for the reads we want to map

    logger.info('target_species_group_df[read_id]: {}'.format(target_species_group_df['read_id']))
    reads = FastqRead.objects.filter(run__flowcell_id=flowcell.id,
                                     read_id__in=target_species_group_df["read_id"])
    # The fasta sequence
    fastqs_list = reads.values_list('read_id', 'fastqreadextra__sequence')

    fastq_input = ''

    for fastq in fastqs_list:

        if fastq[1]:

            temp = '>{}\n{}\n'.format(fastq[0], fastq[1])
            fastq_input += temp

    fastq = fastq_input

    # The command to execute minimap
    map_cmd = '{} -x map-ont -t 4 --secondary=no {} -'.format(minimap2, minimap2_reference_path)
    # Execute minimap2
    out, err = subprocess.Popen(map_cmd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                stderr=subprocess.PIPE).communicate(input=str.encode(fastq))

    if err:
        logger.info("Flowcell id: {} - Minimap error!! {}".format(flowcell.id, err))
    # Get the output of minimap2
    map_output = out.decode()

    if len(map_output) > 0:

        logger.info('>>> Found mappings')
        logger.info('>>> map_output')
        logger.info(map_output)

        logger.info('>>> minimap2 command:')
        logger.info(map_cmd)

    logger.info("Flowcell id: {} - minimap output {} ".format(flowcell.id, map_output))

    plasmid_red_df = pd.DataFrame()
    plasmid_df = species_regions_df[species_regions_df["gff_line_type"] == "plasmid"]
    # If there is output from minimap2 create a dataframe
    if not plasmid_df.empty:
        logger.info("Flowcell id: {} - Mapping reads to plasmids for species {} ".format(flowcell.id, species))

        plasmid_red_series = plasmid_df.apply(plasmid_mapping, args=(species, fastqs_list, flowcell), axis=1)
        plasmid_red_series.dropna(inplace=True)
        logger.info("Flowcell id: {} - plasmid mapping output is {}".format(flowcell.id, plasmid_red_series.head()))
        logger.info("Flowcell id: {} - plasmid mapping output is {}".format(flowcell.id, type(plasmid_red_series)))
        if not plasmid_red_series.empty:
            try:
                for ind, df_plas in plasmid_red_series.iteritems():
                    if type(df_plas) == np.ndarray:
                        df_append = pd.DataFrame(df_plas)
                        df_append.columns = ['read_id', 'query_seq_len', 'query_start', 'query_end', 'rel_strand',
                                             'target_seq_name', 'target_seq_len', 'target_start', 'target_end',
                                             'num_matching_bases', 'num_matching_bases_gaps', 'mapping_qual',
                                             'type_align']
                        logger.info("df_append is {}".format(df_append))
                        plasmid_red_df = plasmid_red_df.append(df_append, sort=True)
                    else:

                        plasmid_red_df = plasmid_red_df.append(df_plas, sort=True)
                logger.info("Flowcell id: {} - plasmid mapping output is {}".format(flowcell.id, plasmid_red_df.head()))
            except AttributeError as e:
                logger.info(e)
                return

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
    logger.info(red_df.keys())
    # If there is output from the plasmid mapping, put it into
    logger.info("Flowcell id: {} - The target df before plasmid mapping results is {}".format(flowcell.id, red_df))
    if not plasmid_red_df.empty and (type(plasmid_red_df) != pd.core.series.Series):
        red_df = red_df.append(plasmid_red_df, sort=True)

    logger.info("Flowcell id: {} - The target df after plasmid mapping results is {}".format(flowcell.id, red_df))

    logger.info(
        "Flowcell id: {} - This many reads mapped evilly on this reference {} for species {}".format(
            flowcell.id,
            red_df.shape[0],
            species)
    )
    # Non red reads, reads that mapped elsewhere on the reference
    # map_df = map_df[~boolean_df]

    logger.info(
        "Flowcell id: {} - This many reads mapped elsewhere on this reference {}".format(flowcell.id,
                                                                                         map_df.shape[0]))
    # get the info on the mapped reads from the original dataframe
    # TODO this is where we sort the barcodes
    update_num_mapped_df = pd.merge(targets_results_df, map_df, how="inner", left_on="read_id", right_on="read_id")
    path = "/home/rory/data/logs/" + flowcell.name + ".csv"
    update_num_mapped_df.to_csv(path, mode="a")
    update_num_mapped_df_non_barcode = update_num_mapped_df.copy(deep=True)

    update_num_mapped_df_non_barcode["barcode_name"] = "All reads"

    update_num_mapped_df_non_barcode["num_mapped"] = update_num_mapped_df_non_barcode.shape[0]

    update_num_mapped_df_non_barcode.drop_duplicates(subset="name", inplace=True)

    update_num_mapped_df_non_barcode.apply(update_mapped_non_dangerous, args=(task,), axis=1)

    # Update the barcoded mapped but none dangerous reads

    update_num_mapped_df_barcoded = update_num_mapped_df

    gb_bc_nm = update_num_mapped_df_barcoded.groupby(["barcode_name", "name"])

    update_num_mapped_df_barcoded.set_index(["barcode_name", "name"], inplace=True)

    update_num_mapped_df_barcoded["num_mapped"] = gb_bc_nm.size()

    update_num_mapped_df_barcoded.reset_index(inplace=True)

    update_num_mapped_df_barcoded.drop_duplicates(subset=["name", "barcode_name"], inplace=True)

    update_num_mapped_df_barcoded.apply(update_mapped_non_dangerous, args=(task,), axis=1)

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


def update_targets_no_mapping(row, task):
    """
    Update the targets to reflect the number of matches, if they don't map
    :param row: The target dataframe row
    :param task: The task database model object
    :return: None
    """
    obj = MappingResult.objects.get(tax_id=row["tax_id"], task=task, barcode_name=row["barcode_name"])
    obj.num_matches = row["num_matches"]
    obj.sum_unique = row["sum_unique"]
    nm = obj.num_mapped
    obj.proportion_of_classified = row["prop_classed"]
    obj.mapped_proportion_of_classified = round(nm / row["num_matches"], 5)
    logger.info(obj.mapped_proportion_of_classified)
    obj.save()


def calculate_num_matches_update(target_df, task, num_matches_targets_barcoded_df, classed_per_bar):
    """
    Calculate the num_matches_update for any target reads
    :param target_df: The dataframe of the target reads
    :param task: The task model object
    :param num_matches_targets_barcoded_df: The number of matches per barcode for each target species in a dataframe
    :param classed_per_bar: The number of reads classified per barcode as a dict
    :return:
    """
    logger.info(target_df)
    barcode_df = target_df.reset_index()
    all_reads_df = target_df.copy(deep=True)
    all_reads_df.reset_index(inplace=True)
    all_reads_df["barcode_name"] = "All reads"
    all_reads_df.drop_duplicates(subset="name", inplace=True)
    if not num_matches_targets_barcoded_df.empty:
        all_reads_df = pd.merge(all_reads_df, num_matches_targets_barcoded_df, how="left", 
                                on=["tax_id", "barcode_name"])
        all_reads_df["num_matches_y"].fillna(0, inplace=True)

        all_reads_df["sum_unique_y"].fillna(0, inplace=True)
        all_reads_df["num_matches"] = all_reads_df["num_matches_x"] + all_reads_df["num_matches_y"]
        all_reads_df["sum_unique"] = all_reads_df["sum_unique_x"] + all_reads_df["sum_unique_y"]

    all_reads_df["prop_classed"] = all_reads_df["num_matches"].div(
        classed_per_bar["All reads"]).mul(100).round(decimals=3)
    logger.info(all_reads_df)
    all_reads_df.apply(update_targets_no_mapping, args=(task,), axis=1)
    # Barcoded now
    gb_bc = barcode_df.groupby(["barcode_name", "name"])

    barcode_df.set_index(["barcode_name", "name"], inplace=True)
    barcode_df["num_matches"] = gb_bc.size()
    barcode_df["sum_unique"] = gb_bc["unique"].sum()
    barcode_df.reset_index(inplace=True)
    if not num_matches_targets_barcoded_df.empty:

        barcode_df = pd.merge(barcode_df, num_matches_targets_barcoded_df, how="left",
                              on=["tax_id", "barcode_name"])

        barcode_df["num_matches_y"].fillna(0, inplace=True)

        barcode_df["sum_unique_y"].fillna(0, inplace=True)

        barcode_df["num_matches"] = barcode_df["num_matches_x"] + barcode_df["num_matches_y"]

        barcode_df["sum_unique"] = barcode_df["sum_unique_x"] + barcode_df["sum_unique_y"]

    barcode_df["total_in_barcode"] = barcode_df["barcode_name"].map(classed_per_bar)

    barcode_df["prop_classed"] = barcode_df["num_matches"].div(
        barcode_df["total_in_barcode"]).mul(100).round(decimals=3)

    barcode_df.drop_duplicates(subset=["barcode_name", "tax_id"], inplace=True)
    logger.info("Barcode Df is {}".format(barcode_df))
    
    barcode_df.apply(update_targets_no_mapping, args=(task,), axis=1)


def map_the_reads(name_df, task, flowcell, num_matches_targets_barcoded_df, targets, class_per_bar):
    """
    Map the reads is called on the dataframe of targets that have been split out, splits them into a group by name, and
    applies the map_all_the_groups function, which returns any reads that map to target areas in a dataframe.
    :param name_df: The targets dataframe
    :param task: The task model object
    :param flowcell: The flowcell model object
    :param num_matches_targets_barcoded_df: The number of matches per barcode for each target species in a dataframe
    :param targets: The targets in this set of targets
    :param class_per_bar: The number of reads classified per barcode as a dict
    :return:
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
                                                               barcode_name="All reads",
                                                               defaults={"red_reads": 0,
                                                                         "num_mapped": 0}
                                                               )
            if created:
                logger.info("Flowcell id: {} - Result Mapping object created for {}".format(flowcell.id, obj.species))
        return

    # Reset the index
    targets_df.reset_index(inplace=True)

    logger.info("Flowcell id: {} - targets_df shape is {}".format(flowcell.id, targets_df.shape))

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
            barcode_name="All reads",
            defaults={"red_reads": 0,
                      "num_mapped": 0,
                      }
        )
        if created:
            logger.info("Flowcell id: {} - Result Mapping object created for {}".format(flowcell.id, obj.species))

    barcodes = targets_df["barcode_name"].unique()

    ncbi = NCBITaxa()

    for barcode in barcodes:
        if barcode not in list(MappingResult.objects.filter(task=task).values_list(
                "barcode_name", flat=True).distinct()):
            for target in targets:
                tax_id = ncbi.get_name_translator([target])[target][0]
                MappingResult(task=task, barcode_name=barcode, species=target, tax_id=tax_id).save()

    del ncbi

    update_num_matches_df = name_df.copy(deep=True)
    if not update_num_matches_df.empty:
        calculate_num_matches_update(update_num_matches_df, task, num_matches_targets_barcoded_df, class_per_bar)

    # Initialise red_df for red alert target reads
    red_reads_df = pd.DataFrame()

    for name, group in gb:
        red_reads_df = red_reads_df.append(map_all_the_groups(group, name, flowcell,
                                                              target_regions_df, targets_df, task,
                                                              num_matches_targets_barcoded_df,
                                                              )
                                           )

    # If none of the reads mapped to inside the groups
    if red_reads_df.empty:
        logger.info("Flowcell id: {} - No reads mapped to dangerous areas! Better luck next time."
                    .format(flowcell.id))
        return

    # results df contains all details on red reads
    results_df = pd.merge(targets_df, red_reads_df, how="inner", left_on="read_id", right_on="read_id")

    results_df_bar = results_df.copy(deep=True)

    mapping_gb_nb = results_df.groupby("name")

    results_df.set_index("name", inplace=True)

    results_df["red_num_matches"] = mapping_gb_nb.size()

    results_df["barcode_name"] = "All reads"

    results_df.reset_index(inplace=True)

    results_df.drop_duplicates(subset="name", inplace=True)

    if "name" in results_df.keys():
        results_df.rename(columns={"name": "species"}, inplace=True)

    logger.info("Flowcell id: {} - The results of the mapping stage pre prev num matches merge "
                "are {}".format(flowcell.id, results_df[["barcode_name", "species", "red_num_matches"]]))

    logger.info("Flowcell id: {} - The results of the mapping stage pre prev num matches merge "
                "are {}".format(flowcell.id, results_df.keys()))

    results_df.apply(update_mapped_red_values, args=(task, flowcell), axis=1)

    # Barcoded red reads

    mapping_gb = results_df_bar.groupby(["name", "barcode_name"])

    results_df_bar.set_index(["name", "barcode_name"], inplace=True)
    # Get how many reads are mapped to each genome dangerously
    # TODO this is where we would combine them so the proportion is correct

    results_df_bar["red_num_matches"] = mapping_gb.size()
    # Get the sum of the unique reads
    results_df_bar.reset_index(inplace=True)

    results_df_bar.rename(columns={"name": "species"}, inplace=True)

    results_df_bar.drop_duplicates(subset=["species", "barcode_name"], inplace=True)

    results_df_bar.apply(update_mapped_red_values, args=(task, flowcell), axis=1)


def create_donut_data_models(row, task):
    """
    Create and return a series of model objects for bulk creation
    :param row: The row of the dataframe
    :param task: The task objects for this analysis
    :return: A model objects to save into the dataframe
    """
    return DonutData(task=task,
                     num_matches=row["num_matches"],
                     sum_unique=row["sum_unique"],
                     name=row["name"],
                     tax_rank=row["tax_rank"],
                     barcode_name=row["barcode_name"])


def update_donut_data_models(row, task):
    """
    Update existing model objectsin the database
    :param row: The row of the dataframe
    :param task: The task object for this analysis
    :return: Nothing
    """
    DonutData.objects.filter(task=task,
                             name=row["name"],
                             barcode_name=row["barcode_name"]).update(num_matches=row["updated_num_matches"],
                                                                      sum_unique=row["updated_sum_unique"])


def calculate_donut_data(df, lineages_df, flowcell, task, tax_rank_filter):
    """
    Calculate the donut data for the donut chart and table
    :param df: The dataframe of centrifuge output, with barcodes as well
    :param lineages_df: The dataframe of lineages for centrifuge output species in this analysis iteration
    :param task: The task object for this analysis
    :param flowcell: The flowcell
    :param tax_rank_filter: A list of the taxonomic ranks in order from superkingdom to species
    :return donut_df: A dataframe containing the
    """
    data_df = pd.merge(df, lineages_df, left_on="tax_id", right_index=True)
    data_df.set_index(tax_rank_filter, inplace=True)

    logger.info('Flowcell id: {} - Calculating donut data'.format(flowcell.id))
    gb_bc = data_df.groupby("barcode_name")
    donut_df = pd.DataFrame()
    for name, group in gb_bc:
        for tax_rank in tax_rank_filter:
            gb = group.groupby(level=tax_rank)
            series = gb["num_matches"].sum()
            temp_df = pd.DataFrame(series)
            temp_df["sum_unique"] = gb["sum_unique"].sum()
            temp_df["tax_rank"] = tax_rank
            temp_df["barcode_name"] = name
            donut_df = donut_df.append(temp_df)

    donut_df = donut_df[donut_df["num_matches"] > 2]

    donut_df.index.name = "name"

    donut_df.reset_index(inplace=True)

    logger.info('Flowcell id: {} - Bulk inserting new species donut data'.format(flowcell.id))

    prev_donut_df = pd.DataFrame(list(DonutData.objects.filter(task=task).values().distinct()))

    if prev_donut_df.empty:
        donut_models_bulk_create = donut_df.apply(create_donut_data_models, args=(task,), axis=1)
        DonutData.objects.bulk_create(list(donut_models_bulk_create.values))
        return
    donut_to_create_df = donut_df[~donut_df["name"].isin(prev_donut_df["name"])]

    donut_models_bulk_create = donut_to_create_df.apply(create_donut_data_models, args=(task,), axis=1)

    DonutData.objects.bulk_create(list(donut_models_bulk_create.values))

    logger.info('Flowcell id: {} - Updating existing species donut data'.format(flowcell.id))

    prev_df_name_barcodes = list(DonutData.objects.filter(task=task).values_list("name", "barcode_name"))

    prev_links_mask = donut_df[['name', 'barcode_name']].agg(tuple, 1).isin(prev_df_name_barcodes)

    donut_data_to_update_df = donut_df[prev_links_mask]

    combined_df = pd.merge(donut_data_to_update_df, prev_donut_df, how="inner", on=["name", "barcode_name"])
    combined_df["updated_num_matches"] = combined_df["num_matches_x"].add(combined_df["num_matches_y"], fill_value=0)
    combined_df["updated_sum_unique"] = combined_df["sum_unique_x"].add(combined_df["sum_unique_y"], fill_value=0)
    combined_df.apply(update_donut_data_models, args=(task,), axis=1)


def run_centrifuge(flowcell_job_id):
    """

    Returns nothing.
    -------
     Run the centrifuge command, and open the subprocess to run it. Keep the process open and pipe the reads into
     it Write the reads incrementally to the processes stdin, capture the stdout and parse it into a useful format.
     Once parsed it's deposited into the database, updating any existing records. The organising key is a uuid,
     stored in self.foreign key.
    """

    task = JobMaster.objects.get(pk=flowcell_job_id)
    task.running = True
    task.save()

    chunk_size = 2000

    flowcell = task.flowcell

    document_number = FastqRead.objects.filter(run__flowcell=flowcell).count()

    logger.info('Flowcell id: {} - Running centrifuge on flowcell {}'.format(flowcell.id, flowcell.name))
    logger.info('Flowcell id: {} - job_master_id {}'.format(flowcell.id, task.id))
    logger.info('Flowcell id: {} - last_read {}'.format(flowcell.id, task.last_read))
    logger.info('Flowcell id: {} - read_count {}'.format(flowcell.id, task.read_count))

    centrifuge_path = get_env_variable("MT_CENTRIFUGE")
    index_path = get_env_variable("MT_CENTRIFUGE_INDEX")
    cmd = "perl " + centrifuge_path + " -f --mm -x " + index_path + " -"

    # Counter for total centOut
    total_centrifuge_output = 0

    # Instance of the Ncbi taxa class, for taxonomic id manipulation
    ncbi = NCBITaxa()

    # Get the task record from the JobMaster table in the database. Used to separate the results inserted,
    # and tell the Javascript and RunMonitor the task is running and when it is complete

    fastqs = FastqRead.objects.filter(run__flowcell=flowcell, id__gt=int(task.last_read)
                                      ).order_by('id')[:chunk_size]

    if fastqs.count() == 0:
        task.complete = True
        logger.info('Flowcell id: {} - Found 0 reads in the database, for this flowcell. Aborting...'.format(
            flowcell.id, chunk_size)
        )
        task.running = False
        task.save()
        return

    if task.read_count + chunk_size > document_number:
        chunk_size = fastqs.count()
        task.save()
        logger.info('Flowcell id: {} - Chunk size is {}, less than 2000 reads in the database'.format(
            flowcell.id, chunk_size)
        )

    logger.info('Flowcell id: {} - number of reads found {}'.format(flowcell.id, fastqs.count()))

    try:
        metadata = Metadata.objects.get(task=task)

    except (Metadata.DoesNotExist, Metadata.MultipleObjectsReturned):

        metadata = Metadata(

            task=task
        )

        metadata.save()

    fastqs_list = fastqs.values_list('read_id', 'barcode_name', 'fastqreadextra__sequence')

    fastqs_data = "".join([str(">read_id=" + c[0] + ",barcode=" + c[1] + "\n" + c[2] + "\n") for c in fastqs_list
                           if c[2] is not None])

    # Write the generated fastq file to stdin, passing it to the command
    logger.info("Flowcell id: {} - Loading index and Centrifuging".format(flowcell.id))

    # Use Popen to run the centrifuge command
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                stderr=subprocess.PIPE).communicate(input=str.encode(fastqs_data))

    if err:
        logger.info("Flowcell id: {} - Centrifuge error!! {}".format(flowcell.id, err))

    # the number of centrifuge lines
    centrifuge_output = out.decode()

    # total number of lines of centrifuge output dealt with
    total_centrifuge_output += centrifuge_output.count("\n") - 1

    logger.info("Flowcell id: {} - number of centrifuge output lines is {}"
                .format(flowcell.id, total_centrifuge_output))

    # output fields is the column headers for the pandas data frame
    output_fields = ["readID", "seqID", "taxID", "numMatches"]

    # create the DataFrame from the output
    df = pd.read_csv(StringIO(centrifuge_output), sep="\t", usecols=output_fields)

    if df.empty:
        logger.info("Flowcell id: {} - No reads found or no centrifuge output."
                    " Check above for error".format(flowcell.id))
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

    barcodes = list(df["barcode_name"].unique())
    barcodes.append("All reads")
    # Get the reads classified in each barcode into a lookup dictionary
    #  so we can work out proportions inside the barcode
    classified_per_barcode = {
        barcode: CentrifugeOutput.objects.filter(task=task, barcode_name=barcode).aggregate(Sum("num_matches"))[
            "num_matches__sum"] for barcode in barcodes}

    if classified_per_barcode["All reads"] is None:
        classified_per_barcode["All reads"] = 0

    for name, group in gb_bc:
        if classified_per_barcode[name] is None:
            classified_per_barcode[name] = 0
        classified_per_barcode[name] += group.shape[0]

    # TODO future mapping target

    # Unique targets in this set of regions
    temp_targets = MappingTarget.objects.filter(target_set="starting_defaults") \
        .values_list("species", flat=True).distinct()

    logger.info("Flowcell id: {} - Targets are {}".format(flowcell.id, temp_targets))

    name_df = df[df["name"].isin(temp_targets)]

    logger.info("Flowcell id: {} - The dataframe target dataframe is {}".format(flowcell.id, name_df.head(n=20)))

    barcode_df = pd.DataFrame()

    # TODO vectorise these bad boys
    logger.info("Flowcell id: {} - The dataframe shape is {}".format(flowcell.id, df.shape))

    unique_barcode = set([c[1] for c in fastqs_list])

    logger.info("Flowcell id: {} - These are the barcodes in this set {}".format(flowcell.id, unique_barcode))

    df2 = df.reset_index()

    df2 = df2.drop_duplicates(keep="first")
    # Make this the results for all reads in the centrifuge output
    df2["barcode_name"] = "All reads"

    # Include these reads in the lookup dict
    classified_per_barcode["All reads"] += df2.shape[0]

    barcode_df = gb_bc.apply(calculate_barcoded_values, barcode_df, classified_per_barcode)

    barcode_df.reset_index(inplace=True)

    barcode_df.rename(columns={"unique": "sum_unique", 0: "num_matches"}, inplace=True)

    barcode_df.set_index("tax_id", inplace=True)

    df_data = CentrifugeOutput.objects.filter(species__in=temp_targets,
                                              task=task,
                                              barcode_name__in=barcodes).values(
        "species", "tax_id", "num_matches", "barcode_name", "sum_unique")

    num_matches_per_target_df = pd.DataFrame(list(df_data))

    logger.info("Flowcell id: {} - The previous number of matches dataframe is {}"
                .format(flowcell.id, num_matches_per_target_df))

    map_the_reads(name_df, task, flowcell, num_matches_per_target_df, temp_targets, classified_per_barcode)

    # delete these columns, no longer needed
    df.drop(columns=["readID", "seqID", "numMatches", "unique", "barcode_name", "read_id"], inplace=True)

    # remove duplicate rows in the data frame,so we only have one entry for each species
    df.drop_duplicates(keep="first", inplace=True)

    # filter out no ranks, like root and cellular organism
    tax_rank_filter = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

    # drop all duplicated lines to keep only one for each entry, creating atomicity
    df.reset_index(inplace=True)

    df = df.drop_duplicates(keep="first")
    # Make this the results for all reads in the centrifuge output
    df["barcode_name"] = "All reads"

    # Include these reads in the lookup dict
    # query the index get all the objects in the database for this analysis at this point
    queryset = CentrifugeOutput.objects.filter(task=task)
    # Get the previous tax ids in this analysis
    prev_df_tax_ids_barcodes = list(queryset.values_list("tax_id", "barcode_name"))

    # ###### BULK CREATE CENTOUTPUT OBJECTS ########
    # Give each row a flowcell object cell, as part of the flowcell_id series
    df["flowcell"] = flowcell

    # same with task from Jobmaster
    df["task"] = task

    df["proportion_of_classified"] = df["num_matches"].div(
        classified_per_barcode["All reads"]).mul(100).round(decimals=3)

    df.set_index("tax_id", inplace=True)
    # Add the barcode dataframe onto the dataframe, so now we have all reads and barcodes
    df = df.append(barcode_df, sort=True)

    df.reset_index(inplace=True)

    lineages_df = insert_new_lineages(ncbi, df, tax_rank_filter, flowcell)

    # df[[tax_rank_filter]].fillna("Unclassified", inplace=True)
    #
    # df[[tax_rank_filter]].replace("nan", "Unclassified", inplace=True, regex=True)
    #
    # df["proportion_of_classified"].fillna(0, inplace=True)

    # Split out the rows that have a tax_id that doesn't already have results
    prev_links_mask = df[['tax_id', 'barcode_name']].agg(tuple, 1).isin(prev_df_tax_ids_barcodes)
    try:
        calculate_donut_data(df, lineages_df, flowcell, task, tax_rank_filter)
    except KeyError as e:
        logger.info("PROBLEMS {}".format(e))
        return

    # subset where there is so these need dootdating
    cent_to_create_df = df[~prev_links_mask]

    cent_to_create_df = pd.merge(cent_to_create_df, lineages_df, how="inner", left_on="tax_id", right_index=True)
    logger.info("Flowcell id: {} - Centrifuge to update dataframe {}".format(flowcell.id, cent_to_create_df.head()))
    cent_to_create_df["task"] = task

    # apply row wise to append model representation objects into a series
    centrifuge_create_series = cent_to_create_df.apply(create_centrifuge_models, args=(classified_per_barcode,), axis=1)

    logger.info("Flowcell id: {} - Bulk creating CentOutput objects".format(flowcell.id))
    # Bulk create the objects
    CentrifugeOutput.objects.bulk_create(list(centrifuge_create_series.values))
    # Get all the results for the barcoded entries, so each species will have multiple entries under different barcodes
    # TODO Split back into one table in database, barcoded results

    cent_to_update_df = df[prev_links_mask]

    prev_tax_ids_array = cent_to_update_df["tax_id"].values

    previous_centrifuge_data_queryset = CentrifugeOutput.objects.filter(task=task,
                                                                        tax_id__in=prev_tax_ids_array).values()
    # If there are results, we need to update existing results database rows and create new ones for non extant rows
    if not cent_to_update_df.empty:
        # Create a pandas dataframe
        previous_df = pd.DataFrame(list(previous_centrifuge_data_queryset))
        logger.info("Flowcell id: {} - Previous CentOutput data found {}".format(flowcell.id, previous_df.head()))
        cent_to_update_df.set_index(["tax_id", "barcode_name"], inplace=True)
        previous_df.set_index(["tax_id", "barcode_name"], inplace=True)
        cent_to_update_df = cent_to_update_df[~cent_to_update_df.index.duplicated(keep="first")]
        previous_df = previous_df[~previous_df.index.duplicated(keep="first")]
        logger.info("Flowcell id: {} - Centrifuge to update dataframe {}".format(flowcell.id, cent_to_update_df.head()))
        cent_to_update_df["num_matches"] = previous_df["num_matches"] + cent_to_update_df["num_matches"]

        cent_to_update_df["sum_unique"] = previous_df["sum_unique"] + cent_to_update_df["sum_unique"]

        cent_to_update_df["proportion_of_classified"] = cent_to_update_df["num_matches"] \
            .div(classified_per_barcode["All reads"]).mul(100).round(decimals=3)

        # ###### Update existing barcode output ######
        # Apply to update the existing entries in the databases
        cent_to_update_df.reset_index(inplace=True)
        cent_to_update_df.apply(update_centrifuge_output_values, args=(task.id,),
                                axis=1)

    # Update the metadata for the headers

    # Get the start time #TODO remove the unused fields
    start_time = metadata.start_time.replace(tzinfo=None)
    end_time = timezone.now().replace(tzinfo=None)

    metadata.finish_time = str(end_time - start_time)
    metadata.save()

    # Update the jobmaster object fields that are relevant
    task = JobMaster.objects.get(pk=task.id)
    task.running = False

    if fastqs.count() > 0:
        task.last_read = fastqs[chunk_size - 1].id
        logger.info("Inseide if")

    task.read_count = task.read_count + fastqs.count()
    task.save()

    logger.info("Flowcell id: {} - New last_read_id is - {}".format(flowcell.id, task.last_read))
    logger.info("Flowcell id: {} - Total CentOut lines - {}".format(flowcell.id, total_centrifuge_output))
    logger.info("Flowcell id: {} - Finished!".format(flowcell.id))
