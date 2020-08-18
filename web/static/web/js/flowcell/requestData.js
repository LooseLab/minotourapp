function requestData (flowcell_id) {
  console.log(`Requesting data`)
  var selected_barcode = getSelectedBarcode(`Metagenomics`)
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

    if (flowcellSelectedTabInput == `summary-data`) {
      // if we are on the summary page - cal l requestRunDetail from the monitor app scope
      // appends the Run summary table HTML onto the summary page
      this.requestRunDetails(flowcell_id)
      /// / append a messages html tableonto the messages div of the summaries tab
      requestMinknowMessages(flowcell_id, data)
    } else if (flowcellSelectedTabInput == `reads`) {
      requestReadData(flowcell_id)
    } else if (flowcellSelectedTabInput == `live-event-data`) {
      // this.requestRunDetails(flowcellId);
      this.requestLiveRunStats(flowcell_id)
    } else if (flowcellSelectedTabInput == `metagenomics`) {
      // The intervals for updating the charts are found in the individual files in the vis-d3 directory
      // SO you are on the Metagenomics tab, congratulations!
      // const addBarcodes = this.addMetaBarcodeTabs.bind(this)
      // // if (selected_barcode == null) {
      // //
      // // }
      // const url = `/api/v1/flowcells/` + flowcell_id + `/metagenomic_barcodes`
      //
      // this.drawSankey(flowcell_id)
      //
      // this.metaHeader(flowcell_id)
      // // Draw the donut rank table
      // this.drawDonutRankTable(flowcell_id)
      // // Draw the donut chart
      // this.drawDonut(flowcell_id)
      // // update the total Reads Table
      // this._getTotalReadsTable(flowcell_id)
      // // Draw the alert mapping targets table;
      // this.update_mapping_table(flowcell_id)
      // // draw the simplified results table
      // this.draw_simple_table(flowcell_id)
      //
      // $.get(url, {}, function (result) {
      //   this.barcodes = result.data.sort()
      //
      //   this.tabs = result.tabs
      //   // draw the sankey
      //   addBarcodes(flowcell_id, this.barcodes, this.tabs)
      //   // update the metadata header
      // })
    } else if (flowcellSelectedTabInput == `advanced-sequence-mapping`) {
      // this.requestAdvancedPafData(flowcellId);
      this.drawReadUntilCharts()
    } else if (flowcellSelectedTabInput == `sequence-assembly`) {
      this.requestGfaData(flowcell_id)
    }
  }.bind(this))
};
