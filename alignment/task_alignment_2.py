from __future__ import absolute_import, unicode_literals

from collections import defaultdict
from itertools import groupby
from pathlib import Path
from typing import List, NamedTuple

import numpy as np
from celery.utils.log import get_task_logger
from django.db.models import F

from alignment.mapper import MAP
from alignment.models import PafSummaryCov
from alignment.utils import (
    get_alignment_result_dir,
    get_or_create_array,
    create_paf_summary_cov,
)
from minknow_data.models import Run
from minotourapp.celery import app
from minotourapp.redis import redis_instance
from reads.models import FastqRead, JobMaster, Barcode
from reference.models import ReferenceInfo

# Set up the logger to write to logging file
logger = get_task_logger(__name__)


@app.task
def run_minimap2_alignment(job_master_id, streamed_reads=None):
    """
    Run a minimap2 based alignment.
    Parameters
    ----------
    job_master_id: int
        The Pk of the JobMaster entry in the database
    streamed_reads: list of dict
        A list of read dictionaries, if reads are being streamed straight into the task
    Returns
    -------
    None

    """
    # Get the JobMaster
    job_master = JobMaster.objects.get(pk=job_master_id)
    if not job_master.reference:
        raise ValueError("Missing Reference file. Please sort out.")
    flowcell = job_master.flowcell
    avg_read_length = int(flowcell.average_read_length)
    if avg_read_length == 0:
        logger.error(
            f"Average read length is zero Defaulting to 450, if this is not the first iteration, "
            f"this is an error."
        )
        avg_read_length = 450
    if not streamed_reads and not isinstance(streamed_reads, list):
        # The chunk size of reads from this flowcell to fetch from the database
        # aim for 50 megabases
        desired_yield = 25 * 1000000
        chunk_size = round(desired_yield / avg_read_length)
        logger.info(f"Fetching reads in chunks of {chunk_size} for alignment.")
        fasta_list = FastqRead.objects.filter(
            flowcell_id=flowcell.id, id__gt=job_master.last_read
        )[:chunk_size].values(
            "read_id",
            "barcode_name",
            "sequence",
            "id",
            "run_id",
            "type__name",
            "barcode_id",
            "is_pass",
            read_type_id=F("type_id"),
        )
        if fasta_list:
            last_read = fasta_list[-1]["id"]
        else:
            logger.info("No fastq found. Skipping this iteration.")
            job_master.running = False
            job_master.save()
            return
    else:
        last_read = job_master.last_read
        fasta_list = streamed_reads
    read_count = len(fasta_list)
    # The chunk size of reads from this flowcell to fetch from the database
    logger.info(
        f"Flowcell id: {flowcell.id} - Running minimap2 on flowcell {flowcell.name}"
    )
    logger.debug(f"Flowcell id: {flowcell.id} - reference {job_master.reference.name}")
    logger.debug(f"Flowcell id: {flowcell.id} - last_read {job_master.last_read}")
    logger.debug(f"Flowcell id: {flowcell.id} - read_count {job_master.read_count}")
    logger.debug(f"Flowcell id: {flowcell.id} - number of reads found {read_count}")
    # If we have pulled back reads, call fasta
    if read_count > 0:
        align_reads_factory(job_master.id, fasta_list, MAP)

    # Update the JobMaster with the metadata after finishing this iteration
    redis_instance.decr("minimaptasks")
    job_master = JobMaster.objects.get(pk=job_master_id)
    job_master.running = False
    job_master.last_read = last_read
    job_master.read_count = job_master.read_count + read_count
    job_master.save()


def align_reads_factory(job_master_id, fasta_list: List, super_function):
    """
    Sends reads to an alignment factory to map reads to the reference.
    Parameters
    ----------
    job_master_id: int
        The PK of the Jobmaster entry.
    fasta_list: list
        DataFrame containing all the sequence data that we need.
    super_function: MappingServer
        Mapping server class for managing mappy

    Returns
    -------
    last_read_pk: int
        The ID of the last read that we dealt with.
    """
    bin_width = 10
    # The JobMaster object
    job_master = JobMaster.objects.get(pk=job_master_id)
    # Unpack metadata about this job
    last_read_pk = int(job_master.last_read)
    reference_pk = int(job_master.reference.id)
    # We're starting fresh boys
    last_read_pk = 0 if last_read_pk is None else last_read_pk
    # The reference for this alignment
    reference_info = ReferenceInfo.objects.get(pk=reference_pk)
    # Get all the chromosomes for this reference
    chromosomes = reference_info.reference_lines.all()
    # A dictionary with all the results for each chromosome
    # For each Chromosome in the chromosomes, set the name of the chromosome as key, and the object as the value
    chromosome_dict = {chromosome.line_name: chromosome for chromosome in chromosomes}
    # add reference if it isn't already in there WARNING fails silently
    if not super_function.valid(reference_info.name):
        super_function.add_reference(
            reference_info.name, reference_info.minimap2_index_file_location
        )
    # send the sequence to the mapper
    results = []
    read_id_2_read_info = {}
    for read_dict in fasta_list:
        read_id_2_read_info[read_dict["read_id"]] = read_dict
        results.extend(
            super_function.map_sequence(reference_info.name, read_dict, True)
        )
    if not results:
        raise FileNotFoundError("No mapping results.")
    run_id = Run.objects.get(pk=fasta_list[0]["run"]).runid
    folder_dir = get_alignment_result_dir(run_id, create=True)
    # Rejected barcode primary key unblocked or sequenced bool, True is unblocked
    get_rejected = {
        barcode.id: barcode.name == "Unblocked"
        for barcode in Barcode.objects.filter(
            run_id=fasta_list[0]["run"], name__in={"Sequenced", "Unblocked"}
        )
    }

    # dictionary of barcode -> contig
    class Contig(NamedTuple):
        """
        named tuple with the path to the contig in it
        """

        path: Path = Path("")
        name: str = ""

        def __str__(self):
            return self.name

    barcode_contig_2_summary_orm = {
        (a_s.barcode_id, a_s.chromosome_name, a_s.read_type.id): a_s
        for a_s in PafSummaryCov.objects.filter(job_master=job_master)
    }
    barcode_2_contig_array = defaultdict(lambda: defaultdict(Contig))
    # sort results on barcode, contig so we can group by them, and minimise the number of opening and closing arrays
    results.sort(key=lambda x: (x[-1], x[5]))

    for key, group in groupby(results, lambda x: (x[-1], x[5], x[6])):
        barcode, contig, contig_length = key
        barcode_orm = Barcode.objects.get(pk=barcode)
        if contig not in barcode_2_contig_array[barcode]:
            barcode_2_contig_array[barcode][contig] = Contig(
                path=get_or_create_array(
                    folder_dir,
                    contig,
                    contig_length=int(contig_length),
                    barcode_name=barcode_orm.name,
                    bin_width=bin_width,
                    create=True
                ),
                name=contig,
            )
        path_2_array = barcode_2_contig_array[barcode][contig].path
        mem_map = np.load(path_2_array, mmap_mode="r+")
        for mapping in group:
            read_id = mapping[0]
            is_rejected = get_rejected[read_id_2_read_info[read_id]["rejected_barcode"]]
            # create a summary for this barcode/contig/read_type
            if (
                barcode,
                contig,
                read_id_2_read_info[read_id]["type"],
            ) not in barcode_contig_2_summary_orm:
                paf_summary = create_paf_summary_cov(
                    job_master,
                    read_id_2_read_info[mapping[0]],
                    reference_info,
                    chromosome_dict[contig],
                )
                barcode_contig_2_summary_orm[
                    (barcode, contig, read_id_2_read_info[read_id]["type"],)
                ] = paf_summary
            else:
                paf_summary = barcode_contig_2_summary_orm[
                    (barcode, contig, read_id_2_read_info[read_id]["type"],)
                ]
            # now to slice and add 1
            mapping_start, mapping_end = (
                int(mapping[7]) // bin_width,
                np.ceil(int(mapping[8]) / bin_width).astype(int),
            )
            mem_map[int(is_rejected), mapping_start: mapping_end + 1] += 1
            paf_summary.total_yield += len(mapping[1])
            # multiple mappings can shaft this
            paf_summary.read_count += 1
        paf_summary.average_read_length = round(
            paf_summary.total_yield / paf_summary.read_count
        )
        paf_summary.coverage = round(
            paf_summary.total_yield / chromosome_dict[contig].chromosome_length
        )
        paf_summary.save()
        mem_map.flush()
