import os
import shutil
from shutil import SameFileError

from Bio import SeqIO
from django.conf import settings
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand, CommandError

from reads.models import Flowcell, Run, Barcode, FastqRead, FastqFile, FastqReadType, FastqReadExtra
from reference.models import ReferenceInfo, ReferenceLine

def parse_fastq_description(description):

    descriptiondict = dict()
    descriptors = description.split(" ")
    del descriptors[0]
    for item in descriptors:
        bits = item.split("=")
        descriptiondict[bits[0]] = bits[1]
    return descriptiondict

class Command(BaseCommand):

    help = 'Add fastqreads from a file'

    def add_arguments(self, parser):
        parser.add_argument('filename', type=str)
        parser.add_argument('user', type=str)
        parser.add_argument('flowcell', type=str)
        parser.add_argument('is_pass', type=str)

    def handle(self, *args, **options):

        print("Processing filename {} for flowcell {}".format(options['filename'], options['flowcell']))

        user = User.objects.get(username=options['user'])

        flowcell, created = Flowcell.objects.get_or_create(
            name=options['flowcell'],
            owner=user
        )

        flowcell.save()

        fastqread_type = FastqReadType.objects.get(name='Template')

        for record in SeqIO.parse(options['filename'], "fastq"):

            description_dict = parse_fastq_description(record.description)

            run, created = Run.objects.get_or_create(

                name=description_dict.get('runid', None),
                flowcell=flowcell,
                owner=user
            )
            run.save()

            fastqfile, created = FastqFile.objects.get_or_create(

                name=os.path.basename(options['filename']),
                run=run,
                runid=run.runid,
                owner=user
            )
            fastqfile.save()

            fastq_read = {}

            fastq_read['read'] = description_dict.get('read', None)
            fastq_read['runid'] = description_dict.get('runid', None)
            fastq_read['channel'] = description_dict.get('ch', None)
            fastq_read['start_time'] = description_dict.get('start_time', None)
            fastq_read['is_pass'] = options['is_pass']
            fastq_read['read_id'] = record.id
            fastq_read['sequence_length'] = len(str(record.seq))

            barcode, created = Barcode.objects.get_or_create(
                name=description_dict.get('barcode', 'No barcode'),
                run=run
            )

            barcode.save()

            # fastqread, created = FastqRead.objects.get_or_create(
            fastqread = FastqRead(

                run=run,
                read_id=fastq_read['read_id'],
                read=fastq_read['read'],
                channel=fastq_read['channel'],
                barcode=barcode,
                barcode_name=barcode.name,
                quality_average=0,
                is_pass=fastq_read['is_pass'],
                type=fastqread_type,
                start_time=fastq_read['start_time']
            )

            fastqread.save()

            fastqread_extra = FastqReadExtra(

                fastqread=fastqread,
                sequence=str(record.seq),
                quality=record.format('fastq').split('\n')[3]
            )

            # fastq_read['fastqfile'] = fastqfile["url"]

            print(record.id, len(record.seq))
