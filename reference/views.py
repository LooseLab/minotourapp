from pathlib import Path

import pyfastx
from django.db.models import Q
from django.db.utils import IntegrityError
from django.shortcuts import render
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from reference.forms import ReferenceForm
from reference.models import ReferenceInfo, ReferenceLine
from reference.serializers import ReferenceInfoSerializer
from reference.utils import validate_reference_checks, create_minimap2_index


@api_view(["GET", "POST", "DELETE"])
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
    if request.method == "GET":
        if request.GET.get("names", False):
            data = ReferenceInfo.objects.filter(
                Q(private=False) | Q(uploader=request.user)
            ).values_list("name", "id")
        else:
            queryset = ReferenceInfo.objects.filter(
                Q(private=False) | Q(uploader=request.user)
            )
            serializer = ReferenceInfoSerializer(
                queryset, many=True, context={"request": request}
            )
            data = serializer.data
            for d in data:
                if d["uploader_id"] == request.user.id:
                    d["deletable"] = True
                else:
                    d["deletable"] = False
        return Response({"data": data}, status=status.HTTP_200_OK)
    elif request.method == "POST":
        form = ReferenceForm(request.POST, request.FILES)
        if form.is_valid():
            files = request.FILES.getlist("file_location")
            for index, file in enumerate(files):
                file_name = Path(file.name)
                name = str(file_name.stem).partition(".")[0]
                # get fastq or fasta
                handle = (
                    pyfastx.Fasta
                    if set(file_name.suffixes).intersection(
                        {".fna", ".fa", ".fsa", ".fasta"}
                    )
                    else pyfastx.Fastq
                )
                # build minimap2 index
                duplicated, sha256_hash = validate_reference_checks(file, request.user)
                if duplicated and len(files) == 1:
                    # ref_info.delete()
                    return Response(f"Exact reference already exists. Please use {sha256_hash}", status=status.HTTP_403_FORBIDDEN)
                elif duplicated:
                    continue
                ref_info = ReferenceInfo(
                    file_name=file.name,
                    name=name,
                    file_location=file,
                    uploader=request.user,
                )
                try:
                    ref_info.save()
                except IntegrityError as e:
                    if len(files) > 1:
                        print(f"File is duplicated! {file.name}")
                    else:
                        return Response(
                            f"Duplicate upload attempted. File - {file.name}",
                            status=status.HTTP_403_FORBIDDEN,
                        )
                minimap2_index_file_location = create_minimap2_index(ref_info, file_name)
                fa = handle(ref_info.file_location.path)
                # Create the Reference info entry in the database
                ref_info, created = ReferenceInfo.objects.update_or_create(
                    file_name=file.name,
                    defaults={
                        "length": fa.size,
                        "minimap2_index_file_location": minimap2_index_file_location,
                        "sha256_checksum": sha256_hash
                    },
                )
                # Create a Reference line entry for each "Chromosome/line"
                for contig in fa:
                    ReferenceLine.objects.create(
                        reference=ref_info,
                        line_name=contig.name,
                        chromosome_length=len(contig),
                    )
            return Response("Hello")
        else:
            return Response(
                "INVALID FORM?", status=status.HTTP_500_INTERNAL_SERVER_ERROR
            )
    elif request.method == "DELETE":
        reference_id = request.data.get("referencePk", None)
        if reference_id:
            ReferenceInfo.objects.get(pk=int(reference_id)).delete()
            return Response("Reference deleted.", status=status.HTTP_200_OK)
        else:
            return Response("Error", status=status.HTTP_503_SERVICE_UNAVAILABLE)
        # Let's delete a reference, with mappings
        pass


@api_view(["GET"])
def reference_manager(request):
    """

    Parameters
    ----------
    request: rest_framework.request.Request
        The request object of the AJAX request
    Returns
    -------

    """
    return render(request, "reference_manager.html", context={"request": request})
