"""
add_gff3_set - add a set of locations from a gff file as target regions for the metagenomics analysis
"""
import pandas as pd
from django.core.management.base import BaseCommand, CommandError
from ete3 import NCBITaxa
from centrifuge.models import MappingTarget
from reference.models import ReferenceInfo
import numpy as np
from rest_framework.authtoken.models import Token
from pathlib import Path


def find_files_of_type(file_or_directory, file_extensions):
    """Return a list of pathlib.Path of files with chosen extensions
    Parameters
    ----------
    file_or_directory : str
        filepath or a directory
    file_extensions : list
        A list of lowercase file extensions including '.' e.g. ['.txt', '.paf']
    Returns
    -------
    list
        If files with extension are found return list of pathlib.Path objects, otherwise return empty list
    """
    file_or_directory = Path(file_or_directory).expanduser()
    if (
        file_or_directory.is_file()
        and "".join(file_or_directory.suffixes).lower() in file_extensions
    ):
        return [file_or_directory]
    elif file_or_directory.is_dir():
        return [
            x
            for x in file_or_directory.iterdir()
            if "".join(x.suffixes).lower() in file_extensions
        ]
    else:
        return []


def gff_create(row, species_name, tax_id, set_name, user_id, private):
    """
    Create the model objects for the gff
    :param row: The row of the data-frame
    :param species_name: The name of the species these gff regions relate to
    :param tax_id: A dictionary with the tax_id as value, species name as key
    :param set_name: The name of the ste to include these reads in
    :param user_id: The id of the user in the database
    :param private: Whether or not to keep this target set private
    :return:
    """
    #
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
        user_id=user_id,
        private=private,
        defaults={"name": name},
    )

    if created:
        print(
            "\033[1;36;1m Adding Entry {} for position {} to {} for species {} in set {}".format(
                name, row["start"], row["end"], species_name, "starting_defaults"
            )
        )
    else:
        print(
            "\033[1;31;1m Entry matching position {} to {} for species {} in set {}"
            " already exists in database".format(
                row["start"], row["end"], species_name, "starting_defaults"
            )
        )


class Command(BaseCommand):
    """
    A command for manage.py to add a gff3 file into the database for a species
    """

    help = (
        "Add a custom Gff3 file to the minoTour database. "
        "Contains the regions for the metagenomics mapping."
        "Please namee the file after the species."
        "It is necessary to have a reference for the species already uploaded,"
        " with the exact name as the gff file name. If not already present, please add one with "
        "python manage.py add_reference."
    )

    def add_arguments(self, parser):
        """
        Add the arguments and flags to the command line
        :param parser:
        :return:
        """
        parser.add_argument(
            "gff",
            help="Path to the input gff files or a directory of gffs"
            ", if a directory is given files with the extensions"
            "'.gff' or '.gff3' will be used. Files can be Gzipped. The species name used is the name of the file."
            " Please seperate with  an underscore",
            nargs="+",
        )
        parser.add_argument(
            "-S",
            "--set",
            type=str,
            help="The name of the target set to include the gff regions "
            "for the species in",
        )

        parser.add_argument(
            "-k",
            "--key",
            type=str,
            required=True,
            help="The api key to connect this target"
            " set with your account. Found in the"
            " profile section of your minotour page,"
            " once logged in.",
        )
        parser.add_argument(
            "-p",
            "--private",
            action="store_true",
            help="Whether or not this target_set will be hidden from other users. Default - false",
        )

    def handle(self, *args, **options):
        """
        Handle the command after it is executed
        :param args:
        :param options: The optional argument flags, if present
        :return:
        """
        try:
            # If a species name has been added
            gff_files = []
            # These should be lowercase and include the '.'
            endings = [".gff", ".gff3", ".gff.gz", ".gff3.gz"]
            for file_or_directory in options["gff"]:
                gff_files.extend(find_files_of_type(file_or_directory, endings))
            # remove none from gff_files
            gff_files = list(filter(None.__ne__, gff_files))

            reference_list = list(
                ReferenceInfo.objects.all().values_list("name", flat=True)
            )

            set_name_list = list(
                MappingTarget.objects.all().values_list("target_set", flat=True)
            )

            user_id = Token.objects.get(key=options["key"]).user_id

            private = options["private"]

            for gff_file in gff_files:
                # Otherwise use the name of the Gff file

                species_name = str(gff_file.stem).partition(".")[0].replace(" ", "_")
                # If we have a set name provided
                if options["set"]:
                    set_name = options["set"].replace(" ", "_")
                # Else raise an issue, as we need one
                else:
                    raise CommandError(
                        "A name is required for this set of targets. Please specify one using -S or --set."
                    )
                # If we don't have a reference for this species, raise a command error
                if species_name not in reference_list:
                    raise CommandError(
                        'No matching reference file with name "%s" found in database. Please upload '
                        "a reference with"
                        " the exact species name" % species_name
                    )
                # If we do have a reference, proceed
                else:
                    print(
                        "\033[1;35;1m Reference found for species {} in database".format(
                            species_name
                        )
                    )
                # Check if the set already exists
                if set_name in set_name_list:
                    print(
                        "\033[1;35;1m Trying to add regions for species {} to already existing set {}".format(
                            species_name, set_name
                        )
                    )
                # Else create a new list
                else:
                    print(
                        "\033[1;35;1m Creating new set {} and adding target regions for species {}".format(
                            set_name, species_name
                        )
                    )

                species_name = species_name.replace("_", " ")

                print(
                    "\033[1;37;1m Processing gff3 file {} with a set name of {}".format(
                        options["gff"], species_name
                    )
                )
                # Instantiate the NCBITaxa for lookup
                ncbi = NCBITaxa()
                # Get the tax_id of this species
                tax_id = ncbi.get_name_translator([species_name])
                # Create a gff dataframe for this file
                gff_df = pd.read_csv(
                    gff_file,
                    sep="\t",
                    header=None,
                    names=[
                        "seq_id",
                        "source",
                        "type",
                        "start",
                        "end",
                        "score",
                        "strand",
                        "phase",
                        "attributes",
                    ],
                    index_col=False,
                )
                # Remove metadata lines
                gff_df = gff_df[~gff_df["seq_id"].str.contains("##")]
                # The name of the target if provided
                gff_df["name"] = (
                    gff_df["attributes"]
                    .str.extract(r"NAME\=(.*)")[0]
                    .str.split(";", expand=True)[0]
                )
                # If no provided name, add this placeholder
                gff_df["name"] = gff_df["name"].fillna("no_provided_name")
                # appl the gff_create function to all rows of the dataframe
                gff_df.apply(
                    gff_create,
                    args=(species_name, tax_id, set_name, user_id, private),
                    axis=1,
                )

        except Exception as e:
            raise CommandError(repr(e))
