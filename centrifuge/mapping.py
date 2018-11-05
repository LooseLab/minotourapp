"""
mapping.py - perform the mappings
"""
from centrifuge.models import CentOutput, CartographyMapped, RedReadIds, CartographyGuide
from reads.models import FastqRead
from reference.models import ReferenceInfo
from django.conf import settings
import pandas as pd
import numpy as np
from io import StringIO
import subprocess
from celery.utils.log import get_task_logger
from reads.models import Flowcell
from jobs.models import JobMaster
pd.options.mode.chained_assignment = None
logger = get_task_logger(__name__)


def map_all_the_groups(group_df, reference_location, flowcell_id):
    """

    :param group_df:
    :param reference_location:
    :param flowcell_id:
    :return:
    """
    species = group_df["name"].unique()
    print("species is {}".format(species))
    species = species[0].replace(" ", "_")
    refs = ReferenceInfo.objects.get(name=species)
    minimap2_ref = reference_location + refs.filename
    reads = FastqRead.objects.filter(run__flowcell_id__in={flowcell_id}, read_id__in=group_df["read_id"])

    sequence_data = list(reads.values_list("fastqreadextra__sequence", flat=True))
    # create list of read ids for zipping
    read_ids = list(reads.values_list("read_id", flat=True))
    # Create list of tuples where
    # the 1st element is the read_id and the second is the sequence
    tupley_list = list(zip(read_ids, sequence_data))
    fastq = "".join([str(">" + c[0] + "\n" + c[1] + "\n")
                     for c in tupley_list if c[1] is not None])
    map_cmd = '{} -x map-ont -t 4 --secondary=no {} -'.format("minimap2", minimap2_ref)
    out, err = subprocess.Popen(map_cmd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE) \
        .communicate(input=str.encode(fastq))
    map_out = out.decode()
    if map_out:
        map_df = pd.read_csv(StringIO(map_out), sep="\t", header=None)
    else:
        print("Flowcell id: {} - No minimap output".format(flowcell_id))
        return

    columns = ["read_id", "query_seq_len", "query_start", "query_end", "rel_strand", "target_seq_name",
               "target_seq_len", "target_start", "target_end", "num_matching_bases", "num_matching_bases_gaps",
               "mapping_qual", "type_align", "number_minimiser", "chaining_score", "chaining_score_2nd_chain",
               "random"]
    map_df.columns = columns
    map_df.drop(columns=["number_minimiser", "chaining_score", "chaining_score_2nd_chain", "random"],
                inplace=True)
    return map_df


def falls_in_region(row, map_df):
    """

    :param row:
    :param map_df:
    :return:
    """
    start_low = map_df["target_start"] > row["start"]
    start_high = map_df["target_start"] < row["end"]
    end_low = map_df["target_end"] > row["start"]
    end_high = map_df["target_end"] < row["end"]

    bool_df = pd.concat([start_low, start_high, end_low, end_high], axis=1, ignore_index=True)
    bool_df["keep"] = np.where((bool_df[0] & bool_df[1]) | (bool_df[2] & bool_df[3]), True, False)
    return bool_df["keep"]


def update_mapped_values(row, bulk_insert_red_id, cart_map_dict, flowcell, task):
    """

    :param row:
    :param bulk_insert_red_id:
    :param cart_map_dict:
    :param flowcell:
    :param task:
    :return:
    """
    if row["name"] not in cart_map_dict:
        mapped = CartographyMapped.objects.get(species=row["name"],
                                               tax_id=row["tax_id"],
                                               flowcell=flowcell, task=task)
        mapped.red_reads = row["red_num_matches"]
        mapped.num_matches = row["num_matches"]
        mapped.sum_unique = row["sum_unique"]

        mapped.save()

        cart_map_dict[row["name"]] = mapped.id

        red_rum = RedReadIds(read_id=row["read_id"], CM_species_id=mapped.id)

        bulk_insert_red_id.append(red_rum)
    else:
        red_rum = RedReadIds(read_id=row["read_id"], CM_species_id=cart_map_dict[row["name"]])

        bulk_insert_red_id.append(red_rum)

    return bulk_insert_red_id


class Metamap:
    """
        Perform the mapping of reads that match our target species
    """
    def __init__(self, dataframe, flowcell_id, task_id, target_set):
        self.dataframe = dataframe
        self.task_id = task_id
        self.target_set = target_set

    def map_the_reads(self):
        """
         Run the map
        :return: The
        """
        # Targets_df
        task = JobMaster.objects.get(pk=self.task_id)
        flowcell = task.flowcell
        targets_df = self.dataframe
        if targets_df.empty:
            return
        targets_df.reset_index(inplace=True)
        # refs = ReferenceInfo.objects.all().values()

        reference_location = getattr(settings, "REFERENCE_LOCATION", None)

        gb = targets_df.groupby(["name"], as_index=False)

        self.target_set = "starting_defaults"

        gff3_df = pd.DataFrame(list(CartographyGuide.objects.filter(set=self.target_set).values()))

        target_species = gff3_df["name"].unique()

        target_tax_id = gff3_df["tax_id"].unqiue()

        target_tuple = list(zip(target_species, target_tax_id))

        for target in target_tuple:
            obj, created = CartographyMapped.objects.get_or_create(flowcell=task.flowcell,
                                                                   task=task,
                                                                   species=target[0],
                                                                   tax_id=target[1],
                                                                   defaults={"red_reads": 0,
                                                                             "num_matches": 0,
                                                                             "sum_unique": 0}
                                                                   )
            if created:
                print("\033[1;36;1m Flowcell id {} - ".format(flowcell.id))
        # All mapped reads
        map_df = gb.apply(map_all_the_groups, reference_location, flowcell.id)

        bool_df = gff3_df.apply(falls_in_region, args=(map_df,), axis=1)
        # Get whether it's inside the boundaries or not
        bool_df = bool_df.any()
        # Red reads
        red_df = map_df[bool_df]

        red_df.set_index(["read_id"], inplace=True)
        # results df contains all details on red reads
        results_df = pd.merge(targets_df, red_df, how="inner", left_on="read_id", right_index=True)
        gb_mp = results_df.groupby(["name"])
        results_df.set_index(["name"], inplace=True)
        results_df["red_num_matches"] = gb_mp.size()
        results_df["sum_unique"] = gb_mp["unique"].sum()
        # TODO this is where barcoding step would be
        results_df.reset_index(inplace=True)
        bulk_insert_red_id = []
        bulk_insert_red_from_update = []
        cart_map_dict = {}
        prev_tax_ids = []
        prev_df = pd.DataFrame(list(CartographyMapped.objects.filter(flowcell__id=flowcell.id,
                                                                     task__id=task.id).values()))

        if not results_df.empty:
            cm_id = prev_df["id"].unique()[0]
            results_df.set_index(["tax_id"], inplace=True)

            results_df["to_add_num_matches"] = prev_df["num_matches"]
            results_df["to_add_red_reads"] = prev_df["red_reads"]
            results_df["to_add_red_sum_unique"] = prev_df["sum_unique"]

            results_df["summed_num_matches"] = results_df["num_matches"] + results_df["to_add_num_matches"]
            results_df["summed_red_reads"] = results_df["red_num_matches"] + results_df["to_add_red_reads"]
            results_df["summed_sum_unique"] = results_df["sum_unique"] + results_df["to_add_red_sum_unique"]

            results_df.reset_index(inplace=True)

            results_df.apply(update_mapped_values, args=(flowcell.id, task.id,
                                                     bulk_insert_red_from_update,
                                                     cm_id), axis=1)
            RedReadIds.objects.bulk_create(bulk_insert_red_from_update)

