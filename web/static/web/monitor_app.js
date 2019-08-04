function set_active_navbar_item(item_index) {

    var nav_bar = document.querySelectorAll('.nav li');

    nav_bar.forEach(function(element, index){

        element.className = '';
        if(index == item_index) element.className += 'active';
    });
}


var FlowcellPageApp = {

    init: function () {
        console.log("initialising flowcellpageapp");
        this.chart_reads_called = null;
        this.chart_yield = null;
        this.chart_average_read_length = null;
        this.chart_maximum_read_length = null;
        //this.average_read_lengths_overtime = null;
        this.average_read_lengths_overtime_new = null;
        this.max_read_lengths_overtime_new = null;
        this.xy_scat_length = null;
        this.trans_top100 = null;
        this.chart_cumulative_number_reads_overtime = null;
        this.chartSequencingRate = null;

        this.barcodes = null;

        this.summaryByMinute = null;
        this.summaryByMinute2 = null;

        this.summary = null;

        this.id = null;
        this.selectedBarcode = null;

        this.rundata = null;

        this.livedatayield = new Array();
        this.livedata = new Array();

        this.livedata.voltage = new Array();
        this.livedata.asictemp = new Array();
        this.livedata.heatsinktemp = new Array();
        this.livedata.strand = new Array();
        this.livedata.adapter = new Array();
        this.livedata.good_single = new Array();
        this.livedata.pore = new Array();
        this.livedata.currpercentage = null;
        this.livedata.currstrand = null;
        this.livedata.percentage = new Array();
        this.livedata.yield_history = new Array();
        this.livedata.meanratio_history = new Array();
        this.livedata.instrand_history = new Array();
        this.livedata.openpore_history = new Array();


        var myStringArray = ["above", "adapter", "below", "good_single", "strand", "inrange", "multiple", "pending_mux_change", "saturated", "unavailable", "unblocking", "unclassified", "unknown", "pore","no_pore","zero"];
        var arrayLength = myStringArray.length;

        this.livedata.pore_history = new Array();

        for (var i = 0; i < arrayLength; i++) {
            this.livedata.pore_history[myStringArray[i]] = new Array();
        };

        this.livedata.minIONname = null;
        this.livedata.colours_string = null;
        this.livedata.scalingfactor = 0;

        this.chart_per_chrom_cov = null;

        this.coveragedata = new Array();
        this.coveragedata.read_type = new Array();

        this.makeColumnChart = makeColumnChart;
        this.makeSplineChart = makeSplineChart;
        this.makeSplineChartNonDatetime = makeSplineChartNonDatetime;
        this.makeBoxPlot = makeBoxPlot;
        this.makeLiveHistogram = makeLiveHistogram;
        this.makeLiveChart = makeLiveChart;
        this.makeAreaPlot = makeAreaPlot;

        this.lastread = 0;

        this.drawSankey = drawSankey;
        this.metaHeader = metaHeader;
        this.drawDonut = drawDonut;
        this.getTotalReadsTable = getTotalReadsTable;
        this.drawDonutRankTable = drawDonutRankTable;
        this.flowcellTaskHistoryTable = flowcellTaskHistoryTable;
        this.addMetaBarcodeTabs = addMetaBarcodeTabs.bind(this);
        this.update_mapping_table = update_mapping_table;
        this.draw_simple_table = draw_simple_table;

        this.drawReadUntilCharts = drawReadUntilCharts;

        this.updatePoreChart = updatePoreChart;

        this.requestChannelSummaryData = requestChannelSummaryData;

        this.requestHistogramData = requestHistogramData.bind(this);

        this.updateAssemblyCharts = updateAssemblyCharts;

        this.updateAssemblyBoxplot = updateAssemblyBoxplot;

        this.createAssemblyTable = createAssemblyTable;

        this.requestGfaData = requestGfaData;

        this.requestPafData = requestPafData.bind(this);

        // this.requestAdvancedPafData = requestAdvancedPafData.bind(this);

        this.requestRunDetails = requestRunDetails.bind(this);

        this.requestLiveRunStats = requestLiveRunStats;

        this.requestSummaryData = requestSummaryData;

        this.requestData = requestData;

        var flowcell_id = get_selected_flowcell();

        this.requestStatistics = requestStatistics.bind(this);

        this.checkFlowcellTabs = checkFlowcellTabs.bind(this);
        this.addStartTabsEvents = addStartTabsEvents.bind(this);
        this.updateBarcodeNavTab = updateBarcodeNavTab.bind(this);

        this.ChartNumContigs = this.makeSplineChartNonDatetime(
            "num-contigs",
            "Number of Contigs Assembled".toUpperCase(),
            "Number of Contigs".toUpperCase(),
            "Number of input Reads".toUpperCase()
        );

        this.ChartN50Contigs = this.makeSplineChartNonDatetime(
            "n50-contigs",
            "Assembly N50".toUpperCase(),
            "Assembly N50".toUpperCase(),
            "Number of input Reads".toUpperCase()
        );

        this.ChartSumContigs = this.makeSplineChartNonDatetime(
            "sum-contigs",
            "Total length of Assembly".toUpperCase(),
            "Total length".toUpperCase(),
            "Number of input Reads".toUpperCase()
        );

        this.ChartBoxPlotContigs = this.makeBoxPlot(
            "contigs-boxplot",
            "Contig Lengths per barcode".toUpperCase(),
            "Contig Length".toUpperCase()
        );

        this.LiveHistogram = this.makeLiveHistogram(
            "live-histogram",
            "Histogram of Read Lengths (Events)",
            "test"
        );

        this.LiveCumuYield = this.makeLiveChart(
            "chart-yield",
            "Yield Over Time",
            "cumulative yield"
        );

        this.LiveInStrand = this.makeLiveChart(
            "live-strand",
            "in strand counts",
            "number of pores in strand/single"
        );

        this.LiveOccupancy = this.makeLiveChart(
            "live-occupancy",
            "% Occupancy Over Time",
            "% Occupancy"
        );

        this.LiveTemperature = this.makeLiveChart(
            "live-temperature",
            "Temperature Over Time",
            "°Celcius"
        );
        this.LiveVoltage = this.makeLiveChart(
            "live-voltage",
            "Voltage Over Time",
            "mV"
        );

        this.PoreShizzle = this.makeAreaPlot(
            "poreshizzle",
            "Pore States".toUpperCase(),
            "Pore States".toUpperCase()
        );


        // this.addStartTabsEvents(flowcell_id);
        // this.checkFlowcellTabs(flowcell_id);

        //('>> calling request data');
        console.log('Calling request data from monitor_app. >>>');
        this.requestData(flowcell_id);
        console.log('Calling request data from monitor_app. <<<');
    }, // end of init

    projectdata: function (data) {
        var results = [];
        var holder = [];
        var diffholder = 0;
        var meanholder = 0;
        if (data.length > 3000) {
            data = data.slice(-3000);
        }


        for (var i = 1; i < data.length; i++) {
            var diff = data[i][1] - data[i - 1][1];
            holder.push(diff);
            meanholder = meanholder + diff;
        }
        for (var i = 2; i < holder.length; i++) {
            var ratio = holder[i] / holder[i - 1];
            //if (ratio > 2){
            //    ratio = 2;
            //}
            diffholder = diffholder + ratio;
            //if (meanholder > 1) {
            //    meanholder = 1;
            //}
        }
        //console.log(diffholder/(holder.length - 1));
        if (diffholder / (holder.length - 1) > 1) {
            return [1, meanholder / (holder.length - 1)];
        } else if (diffholder / (holder.length - 1) < 0.999) {
            return ([0.999, meanholder / (holder.length - 1)]);
        } else {
            return ([diffholder / (holder.length - 1), meanholder / (holder.length - 1)]);
        }
    },

    projectresults: function (syntheticdata, scalingfactor, steps, difference, runstart) {
        testdata = syntheticdata.slice(-2);
        var lastval = testdata[1][1];
        var lasttime = testdata[1][0];
        var timingdiff = testdata[1][0] - testdata[0][0];
        //var valdiff = testdata[1][1]-testdata[0][1];
        var valdiff = difference;
        var newresults = [];
        var muxphase = 0
        //console.log(lastval,timingdiff);
        for (var i = 0; i < steps; i++) {
            templastval = lastval;
            lastval = lastval + (valdiff * scalingfactor);
            valdiff = lastval - templastval;
            lasttime = lasttime + timingdiff;
            if (muxphase != Math.floor((lasttime - (runstart * 1000)) / 1000 / 28800)) {
                difference = difference * 0.9;
                valdiff = difference;
                muxphase = Math.floor((lasttime - (runstart * 1000)) / 1000 / 28800);
            }
            //remainder = lasttime-(runstart*1000) - (lasttime-(runstart*1000) % (28800*1000))/(28800*1000);
            //console.log(Math.floor((lasttime-(runstart*1000))/1000/28800));
            newresults.push([lasttime, Math.ceil(lastval)]);
        }
        //console.log(newresults);
        return newresults;
    },

    scaleyield: function (firstelement, data) {
        var results = [];
        for (var i = 0; i < data.length; i++) {
            //console.log(data[i]);
            //console.log((data[i][0]-firstelement[0])/1000);
            //console.log((data[i][1]/1000000));
            results.push((((data[i][0] - firstelement[0]) / 1000), Math.ceil((data[i][1] / 1000000))));
        }
        return results;
    },


    updateLiveYieldProjection: function () {
        var seqspeed = "450 b/s";
        //commented out below as it wasn't being used - not sure why - need to check
        //var timeleft = this.geteighthours(this.livedatayield.slice(-1), this.rundata.start_time);
        var firsthour = 120;
        [scalingfactor, difference] = this.projectdata(this.livedata.yield_history);
        [scalingfactor2, difference2] = this.projectdata(this.livedata.yield_history.slice(0, firsthour));
        //synthericdata needs renaming!!!
        var syntheticdata = this.livedata.yield_history;
        //console.log(scalingfactor2);
        //console.log(syntheticdata);
        newarray = this.projectresults(syntheticdata, scalingfactor, 4000, difference, new Date(this.rundata.start_time).getTime());
        //console.log(newarray);
        newarray1 = this.projectresults(syntheticdata.slice(0, firsthour), scalingfactor2, 4000, difference2, new Date(this.rundata.start_time).getTime());
        newarray2 = this.projectresults(syntheticdata.slice(0, firsthour), 1, 4000, difference2, new Date(this.rundata.start_time).getTime());
        this.LiveYield.series[0].setData(this.converttobases(this.livedata.yield_history, seqspeed));
        this.LiveYield.series[1].setData(this.converttobases(newarray, seqspeed));
        this.LiveYield.series[2].setData(this.converttobases(newarray1, seqspeed));
        this.LiveYield.series[3].setData(this.converttobases(newarray2, seqspeed));
        this.LiveYield.redraw();
        this.LiveYield.reflow();
    },

    converttobases: function (data, seqspeed) {
        if (Number(this.livedata.scalingfactor) > Number(0)) {
            //console.log("returning scaling factor" + this.livedata.scalingfactor);
            scaling = Number(this.livedata.scalingfactor);
        } else {
            switch (seqspeed) {
                case "MegaCrazy Runs":
                    scaling = 3.5;
                    break;
                case "450 b/s":
                    scaling = 1.8;
                    break;
                case "250 b/s":
                    scaling = 1.1;
                    break;
                case "70 b/s":
                    scaling = 1.0;
                    break;
            }
        }
        var scaleddata = [];
        for (var i = 0; i < data.length; i++) {
            scaleddata.push([data[i][0], data[i][1] * scaling]);
        }
        //console.log("returning estimated scaling factor");
        return scaleddata;
    },
};