import base64
import gzip
import json
import tarfile
from io import BytesIO

import pandas as pd
from celery.utils.log import get_task_logger
from django.shortcuts import render
from django.test import RequestFactory
from matplotlib.figure import Figure
from weasyprint import HTML, CSS

# Create your views here.
from artic.models import ArticBarcodeMetadata
from artic.task_artic_alignment import np
from artic.utils import (
    quick_get_artic_results_directory,
    get_artic_run_stats,
)
from minknow_data.models import Flowcell
from minotourapp.celery import app
from minotourapp.settings import BASE_DIR
from reads.models import JobMaster

logger = get_task_logger(__name__)


def hyphen_to_underscore(dictionary):
    """
    Takes an Array or dictionary and replace all the hyphen('-') in any of its keys with a underscore('_')
    :param dictionary:
    :return: the same object with all hyphens replaced by underscore
    """
    # By default return the same object
    final_dict = dictionary

    # for Array perform this method on every object
    if type(dictionary) is type([]):
        final_dict = []
        for item in dictionary:
            final_dict.append(hyphen_to_underscore(item))

    # for dictionary traverse all the keys and replace hyphen with underscore
    elif type(dictionary) is type({}):
        final_dict = {}
        for k, v in dictionary.items():
            # If there is a sub dictionary or an array perform this method of it recursively
            if type(dictionary[k]) is type({}) or type(dictionary[k]) is type([]):
                value = hyphen_to_underscore(v)
                final_dict[k.replace("-", "_")] = value
            else:
                final_dict[k.replace("-", "_")] = v

    return final_dict


# @app.task
# def build_artic_report(pk: int, request):


@app.task(bind=True)
def build_artic_report(self, pk: int, svg_data):
    """
    Artic report tar ball creation
    Parameters
    ----------
    pk: int
        The primary key of the flowcell the task is on
    svg_data: dict
        The request object of the request that initialised this task - contains svg data
    Returns
    -------
    None

    """
    self.update_state(state="PROGRESS", meta={"done": 100, "total": 200})
    request_factory = RequestFactory()
    my_url = (
        f"api/v1/artic/{pk}/export-report/"  # Replace with your URL -- or use reverse
    )
    request = request_factory.get(my_url)
    flowcell = Flowcell.objects.get(pk=pk)
    task = JobMaster.objects.filter(flowcell_id=pk, job_type_id=16).last()
    tar_file_path = f"/tmp/{task.id}_artic_report.tar.gz"
    tar_file_name = f"{task.id}_artic_report.tar.gz"
    tar_fh = BytesIO()
    with tarfile.open(fileobj=tar_fh, mode="w:gz") as tar:
        barcodes = list(ArticBarcodeMetadata.objects.filter(job_master=task))
        self.update_state(
            state="PROGRESS", meta={"done": 0, "total": len(barcodes) + 3}
        )
        count = 0
        for abm in barcodes:
            count += 1
            data_2 = {}
            selected_barcode = str(abm.barcode.name)
            logger.info(selected_barcode)
            (
                flowcell,
                artic_results_path,
                artic_task_id,
                coverage_path,
            ) = quick_get_artic_results_directory(flowcell.id, selected_barcode)
            try:
                with open(str(coverage_path), "rb") as fh:
                    arr = np.fromfile(fh, dtype=np.uint16)
                    fig = Figure(figsize=(8, 2), dpi=250)
                    ax = fig.subplots()
                    ax.plot(arr)
                    fig.suptitle(f"{selected_barcode} Coverage across genome")
                    ax.set_xlabel(f"Reference position")
                    ax.set_ylabel(f"Coverage")
                    ax.grid()
                    buf = BytesIO()
                    fig.savefig(buf, format="png")
                    data_2[
                        "coverage_graph"
                    ] = f"data:image/png;base64,{base64.b64encode(buf.getvalue()).decode()}"
            except FileNotFoundError as e:
                print(f"Error {repr(e)}")
                self.update_state(
                    state="PROGRESS", meta={"done": count, "total": len(barcodes) + 3}
                )
                continue
            json_path = (
                artic_results_path
                / selected_barcode
                / "json_files"
                / f"{selected_barcode}_ARTIC_medaka.json.gz"
            )
            if json_path.exists():
                with gzip.open(json_path, "r") as fin:
                    data = json.loads(fin.read().decode("utf-8"))
                    data = hyphen_to_underscore(data)
            else:
                data = {}
            csv_path = (
                artic_results_path
                / selected_barcode
                / "csv_files"
                / f"{selected_barcode}_ARTIC_medaka.csv"
            )
            if csv_path.exists():
                df = pd.read_csv(csv_path)
                html_string = df.T.to_html(
                    classes="table table-sm table-responsive", border=0, justify="left"
                )
                data["hidden_html_string"] = html_string

            csv_path = (
                artic_results_path
                / selected_barcode
                / f"{selected_barcode}_ARTIC_medaka.csv.gz"
            )
            if csv_path.exists():
                df = pd.read_csv(csv_path)
                html_string = df.to_html(
                    classes="table table-sm table-responsive",
                    border=0,
                    index=False,
                    justify="left",
                )
                data["hidden_html_string2"] = html_string
            voc_report_html_string = (
                render(
                    request,
                    "artic-variant-of-concern.html",
                    context={"artic_barcode_VoC": data},
                )
                .getvalue()
                .decode()
            )
            data["barcode_name"] = selected_barcode
            data["voc_report_html"] = voc_report_html_string
            # data[
            #     "artic_bar_plot_path"
            # ] = f"file:///{str(artic_results_path)}/{selected_barcode}/{selected_barcode}-barplot.png"
            # data[
            #     "artic_box_plot_path"
            # ] = f"file:///{str(artic_results_path)}/{selected_barcode}/{selected_barcode}-boxplot.png"
            data.update(data_2)
            barcode_report_pdf_bytes = BytesIO(
                HTML(
                    string=render(
                        request, "barcode_report.html", context={"data": data}
                    ).getvalue()
                ).write_pdf(
                    None,
                    stylesheets=[
                        CSS(f"{BASE_DIR}/web/static/web/css/artic-report.css"),
                        CSS(
                            f"{BASE_DIR}/web/static/web/libraries/bootstrap-4.5.0-dist/css/bootstrap.css"
                        ),
                    ],
                )
            )
            barcode_report_pdf_bytes.seek(0, 2)  # go to the end
            source_len = barcode_report_pdf_bytes.tell()
            barcode_report_pdf_bytes.seek(0)
            info = tarfile.TarInfo(f"{task.id}_{selected_barcode}_artic_report.pdf")
            info.size = source_len
            tar.addfile(info, barcode_report_pdf_bytes)
            self.update_state(
                state="PROGRESS", meta={"done": count, "total": len(barcodes) + 3}
            )
        artic_summary_pdf_bytes = get_artic_run_stats(
            pk, svg_data, request, task, logger
        )
        logger.info("Got summary bytes?")
        artic_summary_pdf_bytes.seek(0, 2)  # go to the end
        source_len = artic_summary_pdf_bytes.tell()
        artic_summary_pdf_bytes.seek(0)
        info = tarfile.TarInfo(f"{task.id}_artic_report.pdf")
        info.size = source_len
        tar.addfile(info, artic_summary_pdf_bytes)
        count += 1
        logger.info("trEE SVG")
        self.update_state(
            state="PROGRESS", meta={"done": count, "total": len(barcodes) + 3}
        )
        try:
            tar.add(f"/tmp/{flowcell.id}_tree-plot.pdf", arcname=f"/figures/{flowcell.id}_tree-plot.pdf")
        except FileNotFoundError as e:
            print("tree svg file not found")
        count += 1
        self.update_state(
            state="PROGRESS", meta={"done": count, "total": len(barcodes) + 3}
        )
        logger.info("SNP SVG")
        try:
            print("adding svg")
            tar.add("snp_plot.svg", arcname="/figures/snp_plot.svg")
        except FileNotFoundError as e:
            print("snipit file not found")
        count += 1
        self.update_state(
            state="PROGRESS", meta={"done": count, "total": len(barcodes) + 3}
        )
        logger.info("Writing file")
        with open(f"/tmp/{pk}_artic_report_test.tar.gz", "wb") as f:
            f.write(tar_fh.getvalue())
        logger.info("Done")
        return f"/tmp/{pk}_artic_report_test.tar.gz"
