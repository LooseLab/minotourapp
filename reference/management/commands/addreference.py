import os
import shutil
import subprocess
from shutil import SameFileError

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

            #print(os.path.basename(options['reference']))
            REFERENCELOCATION = getattr(settings, "REFERENCELOCATION", None)

            #print(REFERENCELOCATION)

            srcname = options['reference']
            dstname = os.path.join(REFERENCELOCATION, os.path.basename(options['reference']))

            if srcname.startswith('~') or  dstname.startswith('~'):
                print('Path to reference file and env variable MT_REFERENCE_LOCATION must be absolute.')
                exit()

            #print('origin: {}'.format(srcname))
            #print('destination: {}'.format(dstname))

            try:
                shutil.copy(srcname, dstname)

            except SameFileError as error:
                print('Error trying to copy reference file to Minotour data folder.')
                print(error)
                exit()

            total_length=0
            subfile=dict()
            for record in SeqIO.parse(dstname, "fasta"):
                #print(record.id,len(record.seq))
                subfile[record.id]=len(record.seq)
                total_length=total_length+len(record.seq)
            print("Total Reference Length={}".format(total_length))
            print(subfile)
            cmd = 'bwa index %s' % (dstname)
            subprocess.call(cmd, shell=True)
            cmd2 = "minimap2 -x map-ont -d {}/{}.mmi {}".format(REFERENCELOCATION,os.path.basename(options['reference']),dstname)
            subprocess.call(cmd2, shell=True)
            bwa_index_path = os.path.basename(options['reference'])
            minimap2_index_path = os.path.basename(options['reference']) + ".mmi"
            #print (bwa_index_path)
            print (minimap2_index_path)

            refInfo, created1 = ReferenceInfo.objects.update_or_create(
                reference_name=os.path.basename(options['reference']).split('.')[0],
                filename=os.path.basename(options['reference']),
                #bwa_index_file_location=bwa_index_path,
                minimap2_index_file_location=minimap2_index_path,
                totalrefleN=total_length
            )
            refInfo.save()

            for line in subfile:
                refLine, created2 = ReferenceLine.objects.update_or_create(
                    reference=refInfo,
                    line_name = line,
                    chromosome_length = subfile[line]
                )
                refLine.save()

        except Exception as e:
            raise CommandError(repr(e))
