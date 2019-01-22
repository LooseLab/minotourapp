from __future__ import absolute_import, unicode_literals

from celery import task
from celery.utils.log import get_task_logger
from django.db.models import Count, Max, Min, Sum

from reads.models import FastqRead, Flowcell, RunSummary

logger = get_task_logger(__name__)


@task()
def update_run_summary():
    """
    Update the run summary to be displayed on the flowcell summary tab. Query the mysql client to
    give the stats you need.
    :return:
    """
    logger.info("Updating run summaries")
    # For each active flowcell in all flowcells
    for flowcell in Flowcell.objects.filter(is_active=True):
        # For all the runs under that flowcell
        logger.info("Flowcell id {} - Updating run summaries for flowcell {}".format(flowcell.id, flowcell.id))

        for run in flowcell.runs.all():

            try:
                # Get the run summary for this run
                run_summary_list = RunSummary.objects.filter(run=run)
                # If there is a run summary object the returned queryset will have a length more than 0
                if len(run_summary_list) > 0:
                    # assign the run summary
                    run_summary = run_summary_list[0]
                # If there isn't a run summary object create one
                else:
                    logger.info(
                        "Flowcell id {} - ELSE No objects in database creating new ones".format(flowcell.id)
                    )
                    run_summary = RunSummary.objects.create(

                        run=run,
                        read_count=None,
                        total_read_length=None,
                        max_read_length=None,
                        min_read_length=None,
                        avg_read_length=None,
                        last_read_start_time=None,
                        first_read_start_time=None,
                    )
            # Except No run summary found. One of these code blocks is unnecessary
            except RunSummary.DoesNotExist:
                logger.info(
                    "Flowcell id {} - EXCEPTION No objects in database creating new ones".format(flowcell.id)
                )
                run_summary = RunSummary.objects.create(

                    run=run,
                    read_count=None,
                    total_read_length=None,
                    max_read_length=None,
                    min_read_length=None,
                    avg_read_length=None,
                    last_read_start_time=None,
                    first_read_start_time=None,
                )
            # If the run summary running field is set to False, that means another celery
            # worker isn't currently updating it
            if not run_summary.running:

                run_summary.running = True

                run_summary.save()
                reads = FastqRead.objects.filter(run=run).filter(id__gt=int(run_summary.last_read))[:20000]

                if reads.count() > 0:
                    logger.info(
                        "Flowcell id {} - New reads seen for run summary".format(flowcell.id)
                    )
                    # Do the calculation
                    result = reads.aggregate(
                        Count('sequence_length'),
                        Min('sequence_length'),
                        Max('sequence_length'),
                        Sum('sequence_length'),
                        #Avg('sequence_length'),
                        Min('start_time'),
                        Max('start_time'),
                        Max('id'),
                    )
                    # If there is a previously set read count, so this is not the first iteration
                    if run_summary.read_count != None:
                        logger.info(
                            "Flowcell id {} - Previous results seen for run summary".format(flowcell.id)
                        )
                        # Sum the read counts
                        run_summary.read_count += result['sequence_length__count']
                        # Sum the total lengths
                        run_summary.total_read_length += result['sequence_length__sum']
                        # Get whether the new batch has a read, save it if so or keep old value if not
                        run_summary.max_read_length = max(run_summary.max_read_length, result['sequence_length__max'])
                        # Get whether the new batch a smaller read, save it if so or keep old value if not
                        run_summary.min_read_length = min(run_summary.min_read_length, result['sequence_length__min'])

                        if len(run.RunDetails.all()):
                            run.start_time = run.RunDetails.last().minKNOW_start_time
                            origin = 'Live data'
                            run.save()
                        else:
                            run.start_time = min(run_summary.first_read_start_time,
                                                                result['start_time__min'])
                            run.save()
                        # run_summary.avg_read_length = result['sequence_length__avg']
                        # Get the date of the first read
                        run_summary.first_read_start_time = min(run_summary.first_read_start_time,
                                                                result['start_time__min'])
                        # Get the last read date
                        run_summary.last_read_start_time = max(run_summary.last_read_start_time,
                                                               result['start_time__max'])
                    else:
                        # If this is the first iteration no need to worry about old values, just set the fields
                        # gets turnt up
                        logger.info(
                            "Flowcell id {} - NO Previous results seen for run summary, setting first values"
                            .format(flowcell.id)
                        )
                        run_summary.read_count = result['sequence_length__count']

                        run_summary.total_read_length = result['sequence_length__sum']

                        run_summary.max_read_length = result['sequence_length__max']

                        run_summary.min_read_length = result['sequence_length__min']
                        # run_summary.avg_read_length = result['sequence_length__avg']
                        run_summary.first_read_start_time = result['start_time__min']


                        if len(run.RunDetails.all()):
                            run.start_time = run.RunDetails.last().minKNOW_start_time
                            origin = 'Live data'
                            run.save()
                        else:
                            run.start_time = result['start_time__min']
                            origin = 'Bsecalled Data'
                            run.save()



                        run_summary.last_read_start_time = result['start_time__max']
                    # Get the largest id in your set of reads, set it as a new last read
                    run_summary.last_read = result['id__max']
                    # Get the avg_read_length
                    run_summary.avg_read_length = round(run_summary.total_read_length/run_summary.read_count)

                run_summary.running = False
                # save the results
                run_summary.save()

