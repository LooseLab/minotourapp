

def update_flowcelldetails(flowcell_id):

    flowcell = Flowcell.objects.get(pk=flowcell_id)

    flowcell_details = flowcell.flowcell_details

    # check if flowcell has barcoded run

    flowcellruns = flowcell.flowcelldetails

    flowcell_details.has_barcode = False

    for flowcellrun in flowcellruns.all():
        run = flowcellrun.run
        if run.is_barcoded:
            flowcell_details.has_barcode = True
            break

    flowcell_details.save()
