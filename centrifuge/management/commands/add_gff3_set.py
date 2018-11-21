"""
add_gff3_set - add a set of locations in a gff file
"""
import os
import pandas as pd
from django.core.management.base import BaseCommand, CommandError
from ete3 import NCBITaxa
from centrifuge.models import MappingTarget
from reference.models import ReferenceInfo
import numpy as np


def gff_create(row, set_name, tax_id):
    """
    Create the model objects for the gff
    :param row:
    :param set_name:
    :param tax_id:
    :return:
    """

    if row["name"] != np.nan:
        name = row["name"]
    else:
        name = "danger_zone"

    obj, created = MappingTarget.objects.get_or_create(
        species=set_name,
        tax_id=tax_id[set_name.replace("_", " ")][0],
        target_set="starting_defaults",
        start=row["start"],
        end=row["end"],
        gff_line_type=row["type"],
        defaults={"name": name}
    )

    if created:
        print("\033[1;36;1m Adding Entry {} for position {} to {} for species {} in set {}"
              .format(name, row["start"], row["end"], set_name, "starting_defaults"))
    else:
        print("\033[1;31;1m Entry matching position {} to {} for species {} in set {}"
              " already exists in database"
              .format(row["start"], row["end"], set_name, "starting_defaults"))


class Command(BaseCommand):
    """
    A command for manage.py to add a gff3 file into the database for a species
    """
    help = 'Add a custom Gff3 file to the minoTour database. ' \
           'Contains the target regions for the metagenomics mapping.' \
           'Please either state the species or Name the file after it.'

    def add_arguments(self, parser):
        """
        Add the arguments and flags to the command line
        :param parser:
        :return:
        """
        parser.add_argument('gff', type=str, help="Absolute path to the gff file containing the set of the regions.")
        parser.add_argument('-s', '--species', type=str, help="The name of the species of this genome,"
                                                              " otherwise the filename is used.")

    def handle(self, *args, **options):
        """
        Handle the command after it is executed
        :param args:
        :param options:
        :return:
        """
        try:
            reference_list = list(ReferenceInfo.objects.all().values_list("name", flat=True))

            if options["species"]:
                set_name = options["species"].replace(" ", "_")
            else:
                set_name = os.path.basename(options['gff']).split(".")[0].replace(" ", "_")

            if set_name not in reference_list:
                raise CommandError('No matching reference file with name "%s" found in database. Please upload '
                                   'a reference with'
                                   ' the exact species name' % set_name)

            else:
                print("\033[1;35;1m Reference found for species {} in database".format(set_name))

            set_name = set_name.replace("_", " ")

            print("\033[1;37;1m Processing gff3 file {} with a set name of {}".format(options['gff'], set_name))

            ncbi = NCBITaxa()

            tax_id = ncbi.get_name_translator([set_name])

            srcname = options['gff']

            gff_df = pd.read_csv(srcname, sep="\t", header=None,
                                 names=["seq_id", "source", "type", "start", "end", "score", "strand", "phase",
                                        "attributes"], index_col=False)

            gff_df = gff_df[~gff_df["seq_id"].str.contains("##")]

            gff_df["name"] = gff_df["attributes"].str.extract(r"NAME\=(.*)")[0].str.split(";", expand=True)[0]

            gff_df["name"] = gff_df["name"].fillna("no_provided_name")

            gff_df.apply(gff_create, args=(set_name, tax_id), axis=1)

        except Exception as e:
            raise CommandError(repr(e))


