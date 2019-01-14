import os
import shutil
from shutil import SameFileError

from Bio import SeqIO
from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from reference.models import ReferenceInfo, ReferenceLine


## ToDo: Add ability to handle gz references.

class Command(BaseCommand):

    help = 'Add a custom reference sequence to the minoTour database. ' \
           'Pass the full path to the reference file.'

    def add_arguments(self, parser):
        parser.add_argument('reference', type=str)

    def handle(self, *args, **options):
        try:
            print("Processing Reference {}".format(options['reference']))

            REFERENCE_LOCATION = getattr(settings, "REFERENCE_LOCATION", None)
            MINIMAP2 = getattr(settings, "MINIMAP2", None)

            srcname = options['reference']
            dstname = os.path.join(REFERENCE_LOCATION, os.path.basename(options['reference']))

            if srcname.startswith('~') or dstname.startswith('~'):
                print('Path to reference file and env variable MT_REFERENCE_LOCATION must be absolute.')
                exit()

            try:
                shutil.copy(srcname, dstname)

            except SameFileError as error:
                print('Error trying to copy reference file to Minotour data folder.')
                print(error)
                exit()

            total_length = 0
            subfile=dict()
            for record in SeqIO.parse(dstname, "fasta"):
                #print(record.id,len(record.seq))
                subfile[record.id]=len(record.seq)
                total_length=total_length+len(record.seq)
            print("Total Reference Length={}".format(total_length))
            print(subfile)

            # cmd2 = "{} -x map-ont -d {}/{}.mmi {}".format(
            #     MINIMAP2,
            #     REFERENCE_LOCATION,
            #     os.path.basename(options['reference']),
            #     dstname
            # )

            # subprocess.call(cmd2, shell=True)

            minimap2_index_path = os.path.basename(options['reference']) + ".mmi"

            reference_info, created1 = ReferenceInfo.objects.update_or_create(
                name=os.path.basename(options['reference']).split('.')[0],
                filename=os.path.basename(options['reference']),
                minimap2_index_file_location=minimap2_index_path,
                length=total_length
            )
            reference_info.save()

            for line in subfile:
                reference_line, created2 = ReferenceLine.objects.update_or_create(
                    reference=reference_info,
                    line_name=line,
                    chromosome_length=subfile[line]
                )
                reference_line.save()

        except Exception as e:
            raise CommandError(repr(e))
