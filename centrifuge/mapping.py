from centrifuge.models import CentOutput, CartographyMapped, RedReadIds
from reads.models import FastqRead
from reference.models import ReferenceInfo
from django.conf import settings
import pandas as pd
import numpy as np
from io import StringIO
import subprocess


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
        temp_targets = ["Escherichia coli", "Bacillus cereus"]

        targets_df = self.dataframe

        # refs = ReferenceInfo.objects.all().values()

        reference_location = getattr(settings, "REFERENCELOCATION", None)

        gb = targets_df.groupby(["name"], as_index=False)

        gff3_df = pd.read_csv("/home/rory/data/gff/Escherichia_coli.gff3", sep="\t", header=None,
                              names=["seq_id", "source", "type", "start", "end", "score", "strand", "phase",
                                     "attributes"], index_col=False)

        def map_all_the_groups(group_df):
            species = group_df["name"].unique()
            species = species[0].replace(" ", "_")
            refs = ReferenceInfo.objects.get(reference_name=species)
            minimap2_ref = reference_location + refs.minimap2_index_file_location
            reads = FastqRead.objects.filter(run__flowcell_id__in={self.flowcell_id}, read_id__in=group_df["read_id"])

            sequence_data = list(reads.values_list("fastqreadextra__sequence", flat=True))
            # create list of read ids for zipping
            read_ids = list(reads.values_list("read_id", flat=True))
            # Create list of tuples where
            # the 1st element is the read_id and the second is the sequence
            tupley_list = list(zip(read_ids, sequence_data))
            fastq = "".join([str(">" + c[0] + "\n" + c[1] + "\n")
                             for c in tupley_list if c[1] is not None])
            map_cmd = '{} -x map-ont -t 4 --secondary=no {} -'.format("minimap2", minimap2_ref)
            out, err = subprocess.Popen(map_cmd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)\
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
        map_df = gb.apply(map_all_the_groups)

        def falls_in_region(row):
            start_low = map_df["target_start"] > row["start"]
            start_high = map_df["target_start"] < row["end"]
            end_low = map_df["target_end"] > row["start"]
            end_high = map_df["target_end"] < row["end"]

            bool_df = pd.concat([start_low, start_high, end_low, end_high])
            bool_df["keep"] = np.where((bool_df[0] & bool_df[1]) | (bool_df[2] & bool_df[3]), True, False)
            return bool_df["keep"]
        bool_df = gff3_df.apply(falls_in_region, axis=1)
        bool_df = bool_df.any()
        map_df = map_df[bool_df]
        map_df.set_index(["read_id"])

        results_df = pd.merge(targets_df, map_df, how="inner", left_on="read_id", right_index=True)
        gb_mp = results_df.groupby(["name"])
        results_df["num_matches"] = gb_mp.size()
        results_df["sum_unique"] = gb_mp["unique"].sum()
        # TODO this is where barcoding step would be
        results_df.reset_index(inplace=True)

        results_input_df = results_df.drop_duplicates(subset=["tax_id"])

        bulk_insert = []
        def bulk_create_maps(row):
            mapped = CartographyMapped(species=row["name"], tax_id=row["tax_id"], alert_level=1,
                                       red_reads=row["num_matches"], sum_unique=["sum_unique"])
            bulk_insert.append(mapped)
        # TODO could mark in fastqread? Ask roberto

        def bulk_create_red_reads(row):
            red_read = RedReadIds(read_id=row["read_id"])
















