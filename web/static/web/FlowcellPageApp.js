
var FlowcellPageApp = {

  init: function () {
    console.log(`initialising flowcellpageapp`)
    this.chart_reads_called = null
    this.chart_yield = null
    this.chart_average_read_length = null
    this.chart_maximum_read_length = null
    // this.average_read_lengths_overtime = null;
    this.average_read_lengths_overtime_new = null
    this.max_read_lengths_overtime_new = null
    this.xy_scat_length = null
    this.trans_top100 = null
    this.chart_cumulative_number_reads_overtime = null
    this.chartSequencingRate = null

    this.barcodes = null

    this.summaryByMinute = null
    this.summaryByMinute2 = null

    this.summary = null

    this.id = null
    this.selectedBarcode = null

    this.rundata = null

    this.livedatayield = new Array()
    this.livedata = new Array()

    this.livedata.voltage = new Array()
    this.livedata.asictemp = new Array()
    this.livedata.heatsinktemp = new Array()
    this.livedata.strand = new Array()
    this.livedata.adapter = new Array()
    this.livedata.good_single = new Array()
    this.livedata.pore = new Array()
    this.livedata.currpercentage = null
    this.livedata.currstrand = null
    this.livedata.percentage = new Array()
    this.livedata.yield_history = new Array()
    this.livedata.meanratio_history = new Array()
    this.livedata.instrand_history = new Array()
    this.livedata.openpore_history = new Array()

    var myStringArray = [`above`, `adapter`, `below`, `good_single`, `strand`, `inrange`, `multiple`, `pending_mux_change`, `saturated`, `unavailable`, `unblocking`, `unclassified`, `unknown`, `pore`, `no_pore`, `zero`]
    var arrayLength = myStringArray.length

    this.livedata.pore_history = new Array()

    for (var i = 0; i < arrayLength; i++) {
      this.livedata.pore_history[myStringArray[i]] = new Array()
    };

    this.livedata.minIONname = null
    this.livedata.colours_string = null
    this.livedata.scalingfactor = 0

    this.chart_per_chrom_cov = null

    this.coveragedata = new Array()
    this.coveragedata.read_type = new Array()

    this.makeColumnChart = makeColumnChart
    this.makeSplineChartNonDatetime = makeSplineChartNonDatetime
    this.makeBoxPlot = makeBoxPlot
    this.makeLiveHistogram = makeLiveHistogram
    this.makeLiveChart = makeLiveChart
    this.makeAreaPlot = makeAreaPlot
    this.makePageUnscrollable = makePageUnscrollable

    this.lastread = 0

    this.requestLiveRunStats = requestLiveRunStats

    this.requestData = requestData

    var flowcell_id = getSelectedFlowcell()


    // ('>> calling request data');
    console.log(`Calling request data from monitor_app. >>>`)
    this.requestData(flowcell_id)
    console.log(`Calling request data from monitor_app. <<<`)
  } // end of init
}
