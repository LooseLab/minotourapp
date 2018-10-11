from centrifuge.models import CentOutput, CartographyMapped, RedReadIds
from reads.models import FastqRead
from reference.models import ReferenceInfo
from django.conf import settings
import pandas as pd
import numpy as np
from io import StringIO
import subprocess
from reads.models import Flowcell
from jobs.models import JobMaster


def map_all_the_groups(group_df, reference_location, flowcell_id):
    """

    :param group_df:
    :param reference_location:
    :param flowcell_id:
    :return:
    """
    species = group_df["name"].unique()
    species = species[0].replace(" ", "_")
    refs = ReferenceInfo.objects.get(name=species)
    minimap2_ref = reference_location + refs.minimap2_index_file_location
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
    map_df = pd.read_csv(StringIO(map_out), sep="\t", header=None)

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
    print(bool_df)
    bool_df["keep"] = np.where((bool_df[0] & bool_df[1]) | (bool_df[2] & bool_df[3]), True, False)
    return bool_df["keep"]


def bulk_create_maps(row, bulk_insert_red_id, cart_map_dict, flowcell, task):
    """

    :param row:
    :param bulk_insert_red_id:
    :param cart_map_dict:
    :param flowcell:
    :param task:
    :return:
    """
    if row["name"] not in cart_map_dict:
        mapped = CartographyMapped(species=row["name"], tax_id=row["tax_id"], red_reads=row["red_num_matches"],
                                   num_matches=row["num_matches"], sum_unique=row["sum_unique"],
                                   flowcell=flowcell, task=task)
        mapped.save()
        cart_map_dict[row["name"]] = mapped.id
        red_rum = RedReadIds(read_id=row["read_id"], CM_species_id=mapped.id)
        bulk_insert_red_id.append(red_rum)
    else:
        red_rum = RedReadIds(read_id=row["read_id"], CM_species_id=cart_map_dict[row["name"]])
        bulk_insert_red_id.append(red_rum)

    return bulk_insert_red_id


def to_update_maps(row, flowcell_id, task_id, to_insert_red, cm_id):
    """

    :param row:
    :param flowcell_id:
    :param task_id:
    :param to_insert_red:
    :param cm_id:
    :return:
    """

    CartographyMapped.objects.filter(flowcell__id=flowcell_id, task__id=task_id,
                                     tax_id=row["tax_id"]).update(
        red_reads=row["summed_red_reads"], num_matches=row["summed_num_matches"],
        sum_unique=row["summed_sum_unique"]
    )

    red_id = RedReadIds(read_id=row["read_id"], CM_species_id=cm_id)
    to_insert_red.append(red_id)
    return to_insert_red


class Metamap:
    """
        Perform the mapping of reads that match our target species
    """
    def __init__(self, dataframe, flowcell_id, task_id):
        self.dataframe = dataframe
        self.flowcell_id = flowcell_id
        self.task_id = task_id

    def map_the_reads(self):
        """
         Run the map
        :return: The
        """
        # Targets_df
        flowcell = Flowcell.objects.get(pk=self.flowcell_id)
        task = JobMaster.objects.get(pk=self.task_id)
        targets_df = self.dataframe
        targets_df.reset_index(inplace=True)
        # refs = ReferenceInfo.objects.all().values()

        reference_location = getattr(settings, "REFERENCE_LOCATION", None)

        gb = targets_df.groupby(["name"], as_index=False)

        gff3_df = pd.read_csv("/home/rory/data/gff/Escherichia_coli.gff3", sep="\t", header=None,
                              names=["seq_id", "source", "type", "start", "end", "score", "strand", "phase",
                                     "attributes"], index_col=False)
        # All mapped reads
        map_df = gb.apply(map_all_the_groups, reference_location, self.flowcell_id)

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
        prev_df = pd.DataFrame(list(CartographyMapped.objects.filter(flowcell__id=self.flowcell_id,
                                                                     task__id=self.task_id).values()))
        if not prev_df.empty:
            prev_tax_ids = prev_df["tax_id"]
            prev_df.set_index(["tax_id"], inplace=True)
        to_create_df = results_df[~results_df["tax_id"].isin(prev_tax_ids)]
        to_update_df = results_df[results_df["tax_id"].isin(prev_tax_ids)]

        if not to_create_df.empty:
            to_create_df.apply(bulk_create_maps, args=(bulk_insert_red_id, cart_map_dict,
                                                       flowcell, task), axis=1)
            RedReadIds.objects.bulk_create(bulk_insert_red_id)

        if not to_update_df.empty:
            cm_id = prev_df["id"].unique()[0]
            to_update_df.set_index(["tax_id"], inplace=True)

            to_update_df["to_add_num_matches"] = prev_df["num_matches"]
            to_update_df["to_add_red_reads"] = prev_df["red_reads"]
            to_update_df["to_add_red_sum_unique"] = prev_df["sum_unique"]

            to_update_df["summed_num_matches"] = to_update_df["num_matches"] + to_update_df["to_add_num_matches"]
            to_update_df["summed_red_reads"] = to_update_df["red_num_matches"] + to_update_df["to_add_red_reads"]
            to_update_df["summed_sum_unique"] = to_update_df["sum_unique"] + to_update_df["to_add_red_sum_unique"]

            to_update_df.reset_index(inplace=True)

            to_update_df.apply(to_update_maps, args=(self.flowcell_id, self.task_id,
                                                     bulk_insert_red_from_update,
                                                     cm_id), axis=1)
            RedReadIds.objects.bulk_create(bulk_insert_red_from_update)

