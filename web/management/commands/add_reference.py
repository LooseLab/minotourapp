"""
The django admin management command to add reference files locations to the database,
for mapping and metagenomics validation
"""
import gzip
from pathlib import Path
from Bio import SeqIO
from django.core.management.base import BaseCommand, CommandError
from reference.models import ReferenceInfo, ReferenceLine


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

            for file_or_directory in options["reference"]:
                reference_files.extend(
                    find_files_of_type(file_or_directory, endings)
                )

            # remove none from reference_files
            reference_files = list(filter(None.__ne__, reference_files))
            # TODO doesn't account for private references
            previous_ref = set(
                ReferenceInfo.objects.values_list("name", flat=True).distinct()
            )

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
                    gzip.open(ref_file, "rt")
                    if ref_file.suffix == ".gz"
                    else ref_file
                )
                # Individual lines (I.E Chromosomes in the reference)
                ref_lines = [
                    {"line_name": rec.id, "chromosome_length": len(rec)}
                    for rec in SeqIO.parse(handle, "fasta")
                ]
                # Length of the total reference
                ref_length = sum(
                    [
                        v
                        for r in ref_lines
                        for k, v in r.items()
                        if k == "chromosome_length"
                    ]
                )
                # Create the Reference info entry in the database
                ref_info, created = ReferenceInfo.objects.update_or_create(
                    name=ref_file_stem,
                    filename=ref_file.name,
                    length=ref_length,
                )
                # Create a Reference line entry for each "Chromosome/line"
                for ref_line_dict in ref_lines:
                    ref_line = ReferenceLine(
                        reference=ref_info, **ref_line_dict
                    )
                    ref_line.save()

        except Exception as e:
            raise CommandError(repr(e))
