function requestData (flowcell_id) {
  console.log(`Requesting data`)
  var flowcell_id = flowcell_id

  var url_run = `/api/v1/flowcells/` + flowcell_id + `/`

  $.get(url_run, function (result) {
    // Get data about teflowcell including barcode and number of runs
    const data = result.data

    var barcodes = new Set()
    // For each of the barcodes, add it to the set to be used to draw and select the barcode tabs later
    for (var i = 0; i < data.barcodes.length; i++) {
      if (![`No barcode`, `S`, `U`].includes(data.barcodes[i].name)) {
        barcodes.add(data.barcodes[i].name)
      }
    }

    var flowcellSelectedTabInput = getSelectedTab()
    if (flowcellSelectedTabInput === `reads`) {
      requestReadData(flowcell_id)
    } else if (flowcellSelectedTabInput === `live-event-data`) {
      this.makePageUnscrollable()
      // this.requestRunDetails(flowcellId);
      this.requestLiveRunStats(flowcell_id)
    }
  }.bind(this))
};
