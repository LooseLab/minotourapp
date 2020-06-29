from pathlib import Path

import pyfastx
from django.db.models import Q
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from reference.forms import ReferenceForm
from reference.models import ReferenceInfo, ReferenceLine
from reference.serializers import (ReferenceInfoSerializer)


@api_view(['GET', 'POST'])
def reference_list(request):
    """
    Get a list of all references available to the User making the request.
    Parameters
    ----------
    request: rest_framework.request.Request
        Query params - names. Whether we only want the names and IDs of the references.
    Returns
    -------
    rest_framework.response.Response
    """
    print(request)
    print(request.FILES)
    if request.method == "GET":
        if request.GET.get("names", False):
            data = ReferenceInfo.objects.filter(Q(private=False) | Q(owner=request.user)).values_list("name", "id")
        else:
            queryset = ReferenceInfo.objects.filter(Q(private=False) | Q(owner=request.user))
            serializer = ReferenceInfoSerializer(queryset, many=True, context={'request': request})
            data = serializer.data
        return Response(data, status=status.HTTP_200_OK)
    elif request.method == "POST":
        print(request.data)
        print(request.FILES)
        file_names = request.data.get("filename", [])
        print(file_names)
        if not isinstance(file_names, list):
            file_names = [file_names]
        # md5_checksums = request.data.get("md5", None)
        # files = request.FILES
        # for field_name, file in request.FILES.items():
        #     print(file)
        #     print(dir(file))
        #     print(gzip.open(BytesIO(file.read()), "rt").read())
        #     with open
        form = ReferenceForm(request.POST, request.FILES)
        if form.is_valid():
            files = request.FILES.getlist('file_location')
            print(files)
            for index, file in enumerate(files):
                print(file)
                file_name = Path(file_names[index])
                name = file_name.stem
                ref_info = ReferenceInfo(
                    filename=file_names[0],
                    name=name,
                    file_location=file,
                    owner=request.user
                )
                ref_info.save()
                # If the file is gzipped, unzip it
                handle = (
                    pyfastx.Fastq
                    if file_name.suffix == ".gz"
                    else pyfastx.Fasta
                )
                fa = handle(ref_info.file_location.path)
                # Create the Reference info entry in the database
                ref_info, created = ReferenceInfo.objects.update_or_create(
                    length=fa.size,
                )
                # Create a Reference line entry for each "Chromosome/line"
                for contig in fa:
                    ReferenceLine.objects.create(
                        reference=ref_info, line_name=contig.name, chromosome_length=len(contig)
                    )

            # form.save()
            # print(ReferenceInfo.objects.filter(filename__in=file_names))
            # queryset = ReferenceInfo.objects.filter(filename__in=file_names)
            # for new_ref_file in queryset:
            #     print(new_ref_file.__dict__)
            # # TODO handle file here

            print("HRLLO")
            return Response('Hello')
        else:
            return Response("INVALID FORM?", status=status.HTTP_500_INTERNAL_SERVER_ERROR)