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


def gff_create(row, species_name, tax_id, set_name):
    """
    Create the model objects for the gff
    :param row: The row of the data-frame
    :param species_name: The name of the species these gff regions relate to
    :param tax_id: A dictionary with the tax_id as value, species name as key
    :param set_name: The name of the ste to include these reads in
    :return:
    """

    if row["name"] != np.nan:
        name = row["name"]
    else:
        name = "danger_zone"

    obj, created = MappingTarget.objects.get_or_create(
        species=species_name,
        tax_id=tax_id[species_name.replace("_", " ")][0],
        target_set=set_name,
        start=row["start"],
        end=row["end"],
        gff_line_type=row["type"],
        defaults={"name": name}
    )

    if created:
        print("\033[1;36;1m Adding Entry {} for position {} to {} for species {} in set {}"
              .format(name, row["start"], row["end"], species_name, "starting_defaults"))
    else:
        print("\033[1;31;1m Entry matching position {} to {} for species {} in set {}"
              " already exists in database"
              .format(row["start"], row["end"], species_name, "starting_defaults"))


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
        parser.add_argument("-S", "--set", type=str, help="The name of the target set to include the gff regions "
                                                          "for the species in")
        # TODO add recursive functionality

    def handle(self, *args, **options):
        """
        Handle the command after it is executed
        :param args:
        :param options:
        :return:
        """
        try:
            reference_list = list(ReferenceInfo.objects.all().values_list("name", flat=True))
            set_name_list = list(MappingTarget.objects.all().values_list("target_set", flat=True))

            if options["species"]:
                species_name = options["species"].replace(" ", "_")
            else:
                species_name = os.path.basename(options['gff']).split(".")[0].replace(" ", "_")

            if options["set"]:
                set_name = options["set"].replace(" ", "_")
            else:
                set_name = "starting_defaults"

            if species_name not in reference_list:
                raise CommandError('No matching reference file with name "%s" found in database. Please upload '
                                   'a reference with'
                                   ' the exact species name' % species_name)

            else:
                print("\033[1;35;1m Reference found for species {} in database".format(species_name))

            if set_name in set_name_list:
                print("\033[1;35;1m Trying to add regions for species {} to already existing set {}"
                      .format(species_name, set_name))

            else:
                print("\033[1;35;1m Creating new set {} and adding target regions for species {}"
                      .format(set_name, species_name))

            species_name = species_name.replace("_", " ")

            print("\033[1;37;1m Processing gff3 file {} with a set name of {}".format(options['gff'], species_name))

            ncbi = NCBITaxa()

            tax_id = ncbi.get_name_translator([species_name])

            srcname = options['gff']

            gff_df = pd.read_csv(srcname, sep="\t", header=None,
                                 names=["seq_id", "source", "type", "start", "end", "score", "strand", "phase",
                                        "attributes"], index_col=False)

            gff_df = gff_df[~gff_df["seq_id"].str.contains("##")]

            gff_df["name"] = gff_df["attributes"].str.extract(r"NAME\=(.*)")[0].str.split(";", expand=True)[0]

            gff_df["name"] = gff_df["name"].fillna("no_provided_name")

            gff_df.apply(gff_create, args=(species_name, tax_id, set_name), axis=1)

        except Exception as e:
            raise CommandError(repr(e))
