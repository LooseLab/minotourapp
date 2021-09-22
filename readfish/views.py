from hashlib import sha256

import numpy as np
from django.forms.models import model_to_dict
from django.http import Http404
from django.shortcuts import get_object_or_404
from rest_framework import status
from rest_framework.response import Response
from rest_framework.views import APIView

from artic.models import ArticBarcodeMetadata
from artic.utils import (
    get_amplicon_band_data,
    quick_get_artic_results_directory,
)
from minknow_data.models import Run
from minotourapp.settings import VERSION
from readfish.models import TomlFile, CompletedBarcodes
# Create your views here.
from readfish.utils import get_start_to_start_coords
from reads.models import JobMaster
from reference.models import ReferenceLine


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
            created, t_orm = TomlFile.objects.get_or_create(
                **request.data, sha_hash=self.generate_SHA(request.data["toml_json"])
            )
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

    def get(self, request, run_id, type):
        # run exists
        if type == "run":
            run = get_object_or_404(Run, runid=run_id)
            return Response(model_to_dict(run))
        # job exists
        elif type == "task":
            run = Run.objects.get(runid=run_id)
            jm = get_object_or_404(JobMaster, job_type_id=16, flowcell=run.flowcell)
            return Response(model_to_dict(jm), status=status.HTTP_200_OK)

        else:
            return Response(
                f"Validation of type {type} unknown.",
                status=status.HTTP_400_BAD_REQUEST,
            )


class ConnectHere(APIView):
    """
    Test connection to minoTour for the swordfish client
    """

    def head(self, request):
        return Response(
            status=status.HTTP_200_OK,
            headers={"x-sf-version": ["0.0.1"], "x-mt-ver": VERSION},
        )


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
        CompletedBarcodes.objects.create(run_id=run_id, barcode=barcode_name)

    def get(self, request, run_id, format=None):
        run = Run.objects.get(runid=run_id)
        flowcell = run.flowcell
        job_master = JobMaster.objects.filter(
            job_type_id=16, flowcell_id=flowcell.id
        ).last()
        reference_contig_name = str(
            ReferenceLine.objects.filter(reference_id=job_master.reference_id)
            .first()
            .line_name
        )
        completed_barcodes = set(
            CompletedBarcodes.objects.filter(run_id=run_id).values_list(
                "barcode", flat=True
            )
        )
        barcodes_on_record = set(
            ArticBarcodeMetadata.objects.filter(flowcell=flowcell).exclude(barcode__name="unclassified").values_list(
                "barcode__name", flat=True
            )
        )
        if not barcodes_on_record:
            return Response(
                "No data yet recorded for flowcell", status=status.HTTP_204_NO_CONTENT
            )
        barcodes_to_still_sequence = barcodes_on_record - completed_barcodes
        bed = get_start_to_start_coords(job_master)
        positive_starts = bed[bed["strand"].eq("+")][
            ["range_start", "range_end"]
        ].values.astype(int).tolist()
        negative_starts = bed[bed["strand"].eq("-")][
            ["range_start", "range_end"]
        ].values.astype(int).tolist()
        amplicon_band_coords, _ = get_amplicon_band_data(
            job_master.primer_scheme.scheme_species,
            job_master.primer_scheme.scheme_version,
            job_master.primer_scheme.scheme_directory,
        )
        barcode_rejection_coordinates = {}
        for barcode_name in barcodes_to_still_sequence:
            (
                flowcell,
                artic_results_path,
                artic_task_id,
                coverage_path,
            ) = quick_get_artic_results_directory(flowcell.id, barcode_name)
            try:
                with open(coverage_path, "rb") as fh:
                    coverage = np.fromfile(fh, dtype=np.uint16)
            except FileNotFoundError as e:
                continue
            amplicon_coverages_median = []
            for bin_start, bin_end, _ in amplicon_band_coords:
                amplicon_coverage = coverage[bin_start:bin_end]
                amplicon_median_coverage = np.median(amplicon_coverage)
                amplicon_coverages_median.append(amplicon_median_coverage)
            amplicon_coverages_median = np.array(amplicon_coverages_median)
            amps_to_reject = np.where(amplicon_coverages_median > 20)[0]
            if amps_to_reject.size == len(amplicon_band_coords):
                self.add_barcode_to_completed_table(barcode_name, run_id)
                completed_barcodes.add(barcode_name)
                continue
            positive_barcode_rejection_targets = [
                f"{reference_contig_name},{positive_starts[amp_index][0]},{positive_starts[amp_index][1]},+"
                for amp_index in amps_to_reject
            ]
            negative_barcode_rejection_targets = [
                f"{reference_contig_name},{negative_starts[amp_index][0]},{negative_starts[amp_index][1]},-"
                for amp_index in amps_to_reject
            ]
            positive_barcode_rejection_targets.extend(
                negative_barcode_rejection_targets
            )
            barcode_rejection_coordinates[barcode_name] = {
                "name": barcode_name,
                "control": False,
                "min_chunks": 0,
                "max_chunks": 4,
                "targets": positive_barcode_rejection_targets,
                "single_on": "unblock",
                "single_off": "stop_receiving",
                "multi_on": "unblock",
                "multi_off": "stop_receiving",
                "no_seq": "proceed",
                "no_map": "proceed"
            }
        for completed_barcode_name in completed_barcodes:
            barcode_rejection_coordinates[completed_barcode_name] = {
                "name": completed_barcode_name,
                "control": False,
                "min_chunks": 0,
                "max_chunks": 4,
                "targets": [],
                "single_on": "unblock",
                "single_off": "unblock",
                "multi_on": "unblock",
                "multi_off": "unblock",
                "no_seq": "proceed",
                "no_map": "proceed"
            }
        return Response(barcode_rejection_coordinates, status=status.HTTP_200_OK)
