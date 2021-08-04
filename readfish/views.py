from hashlib import sha256

import numpy as np
from django.http import Http404
from django.shortcuts import get_object_or_404
from rest_framework import status
from rest_framework.response import Response
from rest_framework.views import APIView

from artic.models import ArticBarcodeMetadata
from artic.utils import (
    get_amplicon_band_data,
    slice_coverage_array,
)
from minknow_data.models import Run
from minotourapp.utils import get_env_variable
from readfish.models import TomlFile, CompletedBarcodes
# Create your views here.
from reads.models import JobMaster


class TomlDetail(APIView):
    """
    Retrieve, update or delete a snippet instance.
    """

    def generate_SHA(self, toml_text):
        return sha256(toml_text.encode("utf-8")).hexdigest()

    def get_object(self, run_id):
        try:
            return TomlFile.objects.get(run__runid=run_id)
        except TomlFile.DoesNotExist:
            raise Http404

    def get(self, request, run_id, format=None):
        toml_file_json = self.get_object(run_id).toml_json
        return Response(toml_file_json)

    def post(self, request, run_id, format=None):
        data = request.data
        try:
            created, t_orm = TomlFile.objects.get_or_create(**request.data, sha_hash=self.generate_SHA(request.data["toml_json"]))
            if not created:
                return Response(
                    "Forbidden, already created for this Run Id",
                    status=status.HTTP_403_FORBIDDEN,
                )
        except Exception as e:
            return Response(repr(e), status=status.HTTP_400_BAD_REQUEST)
        return Response(t_orm.__dict__, status=status.HTTP_201_CREATED)


class ValidateSwordfishStartup(APIView):
    """
    Validate whether the task has been created
    """

    def head(self, run_id):
        # run exists
        run = get_object_or_404(Run, runid=run_id)
        # job exists
        jm = get_object_or_404(JobMaster, job_type_id=16, flowcell=run.flowcell)
        return Response(status=status.HTTP_200_OK)


class GetToTheChopper(APIView):
    """
    Put the cookie down
    Return a JSON formatted drop in to the TOML file for updating our targets for readfish barcode
    """

    def add_barcode_to_completed_table(self, barcode_name, run_id):
        """
        Add a barcode to the completed barcodes table
        Parameters
        ----------
        barcode_name: str
            The barcode name

        Returns
        -------
        None
        """
        CompletedBarcodes.objects.create(run_id=run_id, barcode_name=barcode_name)

    def get(self, request, run_id, format=None):
        run = Run.objects.get(runid=run_id)
        flowcell = run.flowcell
        no_overlap_amp_coords, overlap_amp_coords, _ = get_amplicon_band_data(
            get_env_variable("MT_ARTIC_SCHEME_NAME"),
            get_env_variable("MT_ARTIC_SCHEME_VER"),
            no_overlap_only=False,
        )
        completed_barcodes = set(
            CompletedBarcodes.objects.filter(run_id=run_id).values_list(
                "barcode", flat=True
            )
        )
        barcodes_on_record = set(
            ArticBarcodeMetadata.objects.filter(flowcell=flowcell).values_list(
                "barcode__name", flat=True
            )
        )
        if not barcodes_on_record:
            return Response("No data yet recorded for flowcell", status=status.HTTP_204_NO_CONTENT)
        barcodes_to_still_sequence = barcodes_on_record - completed_barcodes
        barcode_coords = {}
        for barcode_name in barcodes_to_still_sequence:
            (
                amplicon_coverages_mean,
                amplicon_coverages_median,
                partial_amplicon_count,
                failed_amplicon_count,
                successful_amplicon_counts,
            ) = slice_coverage_array(no_overlap_amp_coords, flowcell.id, barcode_name)
            median_coverage_array = np.array(amplicon_coverages_median)
            overlap_amp_coords = np.array(overlap_amp_coords)[:, :2]
            targets = overlap_amp_coords[median_coverage_array < 20]
            if targets.size:
                barcode_coords = {barcode_name: targets}
            else:
                self.add_barcode_to_completed_table(barcode_name, run_id)
        return Response(barcode_coords, status=status.HTTP_200_OK)
