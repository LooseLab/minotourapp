"""
The django admin management command to add reference files locations to the database,
for mapping and metagenomics validation
"""
from pathlib import Path

import pyfastx
from django.core.management.base import BaseCommand, CommandError
from rest_framework.authtoken.models import Token

from reference.models import ReferenceInfo, ReferenceLine


def validate_reference_checks(file_path):
    """
    Perform all the checks that we need to do to.
    1. If the reference filename exists in any of the references available to the user, reject it
    2. Check if the md5 of the gzipped file matches any references, if it does, silently backref to that rather
     than add a new file.
    3. If the md5 exists then rename the file to a uuid, move it to the right location and create a md5 checksum for it
    :param file_path: The file path to the reference
    :type file_path: pathlib.Path
    :return:
    """
    pass


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


class Command(BaseCommand):
    """
    The command class from djangos admins interface. Based on argparse
    """

    help = (
        "Add a custom reference sequence to the minoTour database. "
        "Pass the full path to the reference file."
    )
    def add_arguments(self, parser):
        """
        Add the arguments to the base command class
        :param parser:
        :return:
        """
        parser.add_argument(
            "reference",
            help="Path to the input reference files or directories of references"
                 ", if a directory is given files with the extensions"
                 "'.fasta', '.fna' or '.fa' will be used. Files can be gzipped",
            nargs="+",
        )
        parser.add_argument(
            "-k",
            "--key",
            type=str,
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
            Handle the execution of the command
            :param args: The arguments, whether they are present or not
            :param options: The values that have been added to the arguments
            :return:
            """
        try:
            reference_files = []
            # These should be lowercase and include the '.'
            endings = [
                ".fna",
                ".fa",
                ".fasta",
                ".fsa",
                ".fna.gz",
                ".fa.gz",
                ".fasta.gz",
                ".fsa.gz"
            ]
            if options["private"] and not options["key"]:
                print("To add private references, your minotour api_key is required. "
                      "This can be found on the profile page of your account.")
                return
            for file_or_directory in options["reference"]:
                reference_files.extend(
                    find_files_of_type(file_or_directory, endings)
                )
            private = False
            # If we want private references
            if options["private"]:
                private = True
                user = Token.objects.get(key=options["key"]).user
            else:
                user = None
            # remove none from reference_files
            reference_files = list(filter(None.__ne__, reference_files))
            previous_ref = set(
                ReferenceInfo.objects.filter(private=False).values_list("name", flat=True).distinct()
            )
            # If it's private check we aren't multiplying an already existing private reference
            if options["private"]:
                previous_ref = previous_ref.union(set(ReferenceInfo.objects.filter(private=True, owner=user)
                                                      .values_list("name", flat=True).distinct()))
            for ref_file in reference_files:
                # Get the species name of this reference, no file suffixes
                ref_file_stem = str(ref_file.stem).partition(".")[0]
                print("Processing file: {}".format(ref_file.name))
                if ref_file_stem in previous_ref:
                    print(
                        "A reference already exists for this species name: {}".format(
                            ref_file_stem
                        )
                    )
                    print(
                        "If you believe this is in error, or want to add this reference anyway,"
                        " please change the filename"
                    )
                    continue
                # If the file is gzipped, unzip it
                handle = (
                    pyfastx.Fastq
                    if ref_file.suffix == ".gz"
                    else pyfastx.Fasta
                )
                # Individual lines (I.E Chromosomes in the reference)
                fa = handle(ref_file, build_index=False)
                # Create the Reference info entry in the database
                ref_info, created = ReferenceInfo.objects.update_or_create(
                    name=ref_file_stem,
                    filename=ref_file.name,
                    length=fa.size,
                    private=private,
                    owner=user
                )
                # Create a Reference line entry for each "Chromosome/line"
                for contig in fa:
                    ReferenceLine.objects.create(
                        reference=ref_info, line_name=contig.name, chromosome_length=len(contig)
                    )
                print("Successfully handled file.")

        except Exception as e:
            raise CommandError(repr(e))