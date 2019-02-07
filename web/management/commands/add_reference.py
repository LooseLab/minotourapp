import gzip
import os

from Bio import SeqIO
from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from reference.models import ReferenceInfo, ReferenceLine


class Command(BaseCommand):

    help = 'Add a custom reference sequence to the minoTour database. ' \
           'Pass the full path to the reference file.'

    def add_arguments(self, parser):
        parser.add_argument('reference', type=str)

    def handle(self, *args, **options):
        try:
            print("Processing Reference {}".format(options['reference']))

            srcname = options['reference']

            if srcname.startswith('~') or srcname.startswith('~'):
                print('Path to reference file and env variable MT_REFERENCE_LOCATION must be absolute.')
                exit()

            minimap2_index_path = os.path.basename(options['reference']) + ".mmi"

            reference_info, created1 = ReferenceInfo.objects.update_or_create(
                name=os.path.basename(options['reference']).split('.')[0],
                filename=os.path.basename(options['reference']),
                minimap2_index_file_location=minimap2_index_path,
                length=0
            )

            reference_info.save()

            total_length = 0

            if srcname.endswith(".gz"):

                handle = gzip.open(srcname, "rt")

            else:

                handle = srcname

            for record in SeqIO.parse(handle, "fasta"):

                total_length = total_length+len(record.seq)

                reference_line, created2 = ReferenceLine.objects.update_or_create(
                    reference=reference_info,
                    line_name=record.id,
                    chromosome_length=len(record.seq)
                )
                reference_line.save()

            reference_info.length = total_length
            reference_info.save()

            print("Total Reference Length={}".format(total_length))

        except Exception as e:
            raise CommandError(repr(e))

