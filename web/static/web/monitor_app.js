function set_active_navbar_item(item_index) {

    var nav_bar = document.querySelectorAll('.nav li')

    nav_bar.forEach(function(element, index){

        element.className = '';
        if(index == item_index) element.className += 'active';
    });
}


var FlowcellPageApp = {

    constructor () {

    },

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
        this.livedata.good_single = new Array();
        this.livedata.currpercentage = null;
        this.livedata.currstrand = null;
        this.livedata.percentage = new Array();
        this.livedata.yield_history = new Array();
        this.livedata.meanratio_history = new Array();
        this.livedata.instrand_history = new Array();
        this.livedata.openpore_history = new Array();

        var myStringArray = ["above", "adapter", "below", "good_single", "strand", "inrange", "multiple", "pending_mux_change", "saturated", "unavailable", "unblocking", "unclassified", "unknown"];
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

        this.makeChart = makeChart;
        this.makeChart2 = makeChart2;
        this.makeChart3 = makeChart3;
        this.makeChart4 = makeChart4;
        this.makeBoxPlot = makeBoxPlot;
        this.makeChartlabels = makeChartlabels;
        this.makeLiveHistogram = makeLiveHistogram;
        this.makeYieldProjection = makeYieldProjection;
        this.makeLiveChart = makeLiveChart;
        this.makeAreaPlot = makeAreaPlot;

        this.lastread = 0;
        this.needtoupdatecharts = false;

        this.drawSankey = drawSankey;
        this.metaHeader = metaHeader;
        this.drawDonut = drawDonut;
        this.getTotalReadsTable = getTotalReadsTable;
        this.drawDonutRankTable = drawDonutRankTable;
        this.flowcellTaskHistoryTable = flowcellTaskHistoryTable;
        this.addMetaBarcodeTabs = addMetaBarcodeTabs.bind(this);
        this.update_mapping_table = update_mapping_table;

        this.updatePoreChart = updatePoreChart;

        this.requestChannelSummaryData = requestChannelSummaryData;

        this.requestHistogramData = requestHistogramData.bind(this);

        this.updateAssemblyCharts = updateAssemblyCharts;

        this.updateAssemblyBoxplot = updateAssemblyBoxplot;

        this.createAssemblyTable = createAssemblyTable;

        this.requestGfaData = requestGfaData;

        this.requestPafData = requestPafData.bind(this);

        // this.liveUpdateTasks = liveUpdateTasks;
        //
        // this.updateTasks = updateTasks.bind(this);

        // this.requestTasks = requestTasks.bind(this);

        this.requestRunDetails = requestRunDetails.bind(this);

        this.requestLiveRunStats = requestLiveRunStats;

        this.requestSummaryData = requestSummaryData;

        this.requestData = requestData;
        var flowcell_id = get_selected_flowcell();

        this.requestStatistics = requestStatistics.bind(this);

        this.checkFlowcellTabs = checkFlowcellTabs.bind(this);
        this.addStartTabsEvents = addStartTabsEvents.bind(this);
        this.updateBarcodeNavTab = updateBarcodeNavTab.bind(this);

        this.ChartNumContigs = this.makeChart4(
            "num-contigs",
            "Number of Contigs Assembled".toUpperCase(),
            "Number of Contigs".toUpperCase(),
            "Number of input Reads".toUpperCase()
        );

        this.ChartN50Contigs = this.makeChart4(
            "n50-contigs",
            "Assembly N50".toUpperCase(),
            "Assembly N50".toUpperCase(),
            "Number of input Reads".toUpperCase()
        );

        this.ChartSumContigs = this.makeChart4(
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
        /*
        this.LivePoreState = this.makeLiveChart(
            "live-porestate",
            "Pore State Currents",
            "Current pA"
        );
        */

        /*
        this.LiveCurrentRatio = this.makeLiveChart(
            "live-currentratio",
            "Current Ratio In Strand/Open Pore",
            "Current Ratio"
        );
        */

        this.PoreShizzle = this.makeAreaPlot(
            "poreshizzle",
            "Pore States".toUpperCase(),
            "Pore States".toUpperCase()
        );

        // var inputFlowcellId = document.querySelector("#flowcell-id");
        //
        // this.flowcellId = inputFlowcellId.value;

        this.addStartTabsEvents(flowcell_id);
        this.checkFlowcellTabs(flowcell_id);

        console.log('>> calling request data');
        this.requestData(flowcell_id);
    }, // end of init

    updatePoreStats: function () {
        var returndata = this.parseporehist(this.livedata.colours_string, this.livedata.pore_history);
        //console.log(returndata);
        //returndata.sort(a,b){

        //}
        while (this.PoreShizzle.series.length > 0)
            this.PoreShizzle.series[0].remove(true);
        //this.PoreShizzle.addSeries(returndata[4]);
        //console.log(returndata[4]);
        for (var i = 0; i < returndata.length; i++) {
            console.log(returndata[i]);
            var seriesdata = returndata[i];
            //seriesdata['data'].sort(function (a, b){
            //    return a[0] - b[0];
            //});
            this.PoreShizzle.addSeries(seriesdata);
        }
    },

    parseporehist: function (descriptions, counts) {
        var results = [];
        var colors = [];
        var categories = [];
        var datam = [];
        var colorlookup = [];
        //console.log(descriptions);
        descriptions = JSON.parse(descriptions);
        for (var thing in descriptions["groups"]) {
            //console.log(thing);
            if (descriptions["groups"].hasOwnProperty(thing)) {
                if (descriptions["groups"][thing].hasOwnProperty("style")) {
                    //console.log(descriptions["groups"][thing]["states"]);
                    for (var state in descriptions["groups"][thing]["states"]){
                        //console.log(descriptions["groups"][thing]["states"][state]["name"]);
                        colorlookup[descriptions["groups"][thing]["states"][state]["name"]]=descriptions["groups"][thing]["states"][state]["style"]["colour"];
                    }
                    //colorlookup[descriptions["groups"][thing]["name"]] = descriptions["groups"][thing]["style"]["colour"];
                }
            }
        }
        for (var pore in counts) {
            results.push({"name": pore, "color": "#" + colorlookup[pore], "data": counts[pore]})//,"color":"#121212"]});
        }
        return results
    },

    updateLiveCumuYield: function () {
        //console.log(this.LiveCumuYield);
        if (this.LiveCumuYield.series.length < 1) {
            this.LiveCumuYield.addSeries({
                data: this.livedata.yield_history
            });
            this.LiveCumuYield.series[0].update({name: "Events"}, false);
        } else {
            this.LiveCumuYield.series[0].setData(this.livedata.yield_history);
            this.LiveCumuYield.series[0].update({name: "Events"}, false);
        }
        this.LiveCumuYield.redraw();
        this.LiveCumuYield.reflow();
        if (this.LiveOccupancy.series.length < 1) {
            this.LiveOccupancy.addSeries({
                data: this.livedata.percentage
            });
            this.LiveOccupancy.series[0].update({name: "% Occupancy"}, false);
        } else {
            this.LiveOccupancy.series[0].setData(this.livedata.percentage);
            this.LiveOccupancy.series[0].update({name: "% Occupancy"}, false);
        }
        this.LiveOccupancy.redraw();
        this.LiveOccupancy.reflow();
        if (this.LiveInStrand.series.length < 2) {
            this.LiveInStrand.addSeries({data: this.livedata.strand});
            this.LiveInStrand.series[0].update({name: "In Strand"}, false);
            this.LiveInStrand.addSeries({data: this.livedata.good_single});
            this.LiveInStrand.series[1].update({name: "Single Pore"}, false);

        } else {
            this.LiveInStrand.series[0].setData(this.livedata.strand);
            this.LiveInStrand.series[0].update({name: "In Strand"}, false);
            this.LiveInStrand.series[1].setData(this.livedata.good_single);
            this.LiveInStrand.series[1].update({name: "Single Pore"}, false);
        }
        this.LiveInStrand.redraw();
        this.LiveInStrand.reflow();

        /*
        if (this.LivePoreState.series.length < 1) {
            this.LivePoreState.addSeries({data: this.livedata.instrand_history});
            this.LivePoreState.addSeries({data: this.livedata.openpore_history});
        } else {
            this.LivePoreState.series[0].setData(this.livedata.instrand_history);
            this.LivePoreState.series[1].setData(this.livedata.openpore_history);
        }
        this.LivePoreState.series[0].update({name: "In Strand"}, false);
        this.LivePoreState.series[1].update({name: "Open Pore"}, false);
        this.LivePoreState.redraw();
        this.LivePoreState.reflow();
        */

        /*
        if (this.LiveCurrentRatio.series.length < 1) {
            this.LiveCurrentRatio.addSeries({data: this.livedata.meanratio_history});
        } else {
            this.LiveCurrentRatio.series[0].setData(this.livedata.meanratio_history);
        }
        this.LiveCurrentRatio.series[0].update({name: "Current Ratio"}, false);
        this.LiveCurrentRatio.redraw();
        this.LiveCurrentRatio.reflow();
        */
        
        if (this.LiveTemperature.series.length < 1) {
            this.LiveTemperature.addSeries({data: this.livedata.asictemp});
            this.LiveTemperature.addSeries({data: this.livedata.heatsinktemp});
        } else {
            this.LiveTemperature.series[0].setData(this.livedata.asictemp);
            this.LiveTemperature.series[1].setData(this.livedata.heatsinktemp);
        }
        this.LiveTemperature.series[0].update({name: "Asic Temp"}, false);
        this.LiveTemperature.series[1].update({name: "HeatSink Temp"}, false);
        this.LiveTemperature.redraw();
        this.LiveTemperature.reflow();
        if (this.LiveVoltage.series.length < 1) {
            this.LiveVoltage.addSeries({data: this.livedata.voltage});
        } else {
            this.LiveVoltage.series[0].setData(this.livedata.voltage);
        }
        this.LiveVoltage.series[0].update({name: "Voltage"}, false);
        this.LiveVoltage.redraw();
        this.LiveVoltage.reflow();

    },


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
        this.LiveYield.reflow()
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

    tohistogram: function (readeventcountweightedhist, readeventcountweightedhistbinwidth, totalyield) {
        var results = [];
        var categories = [];
        //var counter = 0;
        readeventcountweightedhist = readeventcountweightedhist.replace(/u(?=[^:]+')/g, "").replace(/'/g, "");
        readeventcountweightedhist = JSON.parse(readeventcountweightedhist);
        //console.log(readeventcountweightedhist);
        var n50count = 0;
        var n50index = 0;
        var check = 0;
        //console.log(readeventcountweightedhist);
        for (i in readeventcountweightedhist) {
            //if (readeventcountweightedhist[i] > 0){
            //counter+=1;
            //console.log(readeventcountweightedhistbinwidth);
            //console.log(i);

            //console.log(i*readeventcountweightedhistbinwidth, readeventcountweightedhist[i]);
            n50count += parseInt(readeventcountweightedhist[i]);
            if (n50count >= (parseInt(totalyield) / 2)) {
                //console.log("n50",(i+1)*readeventcountweightedhistbinwidth, n50count);
                check += 1;
            }
            //console.log(i);
            //console.log(parseInt(i)+1);
            var category = String((parseInt(i)) * readeventcountweightedhistbinwidth) + " - " + String((parseInt(i) + 1) * readeventcountweightedhistbinwidth) + " ev";
            categories.push(category);
            if (check == 1) {
                n50index = i;
                results.push({"name": category, "y": parseInt(readeventcountweightedhist[i]), "color": "red"});
                check += 1;
            } else {
                results.push({"name": category, "y": parseInt(readeventcountweightedhist[i]), "color": "blue"});
            }

            //}
        }
        categories.push(">> max ev");
        var missed = 0;
        //var missed = totalyield - readeventcountweightedhist.reduce(add,0);
        results.push({"name": ">> max ev", "y": missed});
        //console.log(n50index);
        return [results, categories, n50index];

    },

    updateLiveHistogram: function (data) {
        returndata = this.tohistogram(data[data.length - 1].minKNOW_histogram_values, data[data.length - 1].minKNOW_histogram_bin_width, data[data.length - 1].event_yield);
        this.LiveHistogram.series[0].setData(returndata[0]);
        this.LiveHistogram.xAxis[0].setCategories(returndata[1]);
        var N50 = parseInt(returndata[2]);
        this.LiveHistogram.xAxis[0].removePlotBand("plot-band-1");
        this.LiveHistogram.xAxis[0].addPlotBand({
            from: N50 - 0.5,
            to: N50 + 0.5,
            color: "#FCFFC5",
            id: "plot-band-1",
        });
        this.LiveHistogram.xAxis[0].removePlotBand("plot-band-2");
        this.LiveHistogram.xAxis[0].addPlotBand({
            color: "black",
            width: 2,
            dashStyle: "longdashdot",
            value: returndata[2],
            label: {
                text: "Estimated Read N50",
                align: "left",
                rotation: 0,
                x: +10 // Amount of pixels the label will be repositioned according to the alignment.
            },
            id: "plot-band-2",
        });
        this.LiveHistogram.xAxis[0].removePlotBand("plot-band-3");
        this.LiveHistogram.xAxis[0].addPlotBand({
            color: "black",
            width: 2,
            dashStyle: "longdashdot",
            value: (Math.floor(this.totalyield / this.readcount / this.datain2)),
            label: {
                text: "Estimated Read Average - " + Math.round(this.totalyield / this.readcount / 1000 * 100) / 100 + " K events",
                align: "left",
                rotation: 0,
                x: +10,
                y: +30, // Amount of pixels the label will be repositioned according to the alignment.
            },
            id: "plot-band-3",
        });
        this.LiveHistogram.reflow();
    },


    calculatereadtoeventscaling: function () {
        var totalyield = 0;
        var readcount = 0;
        if (this.summary !== null) {
            //if ("All reads" in this.summary) {
            for (var readtype in this.summary["All reads"]) {
                //console.log(this.summary["All reads"][readtype]);
                totalyield = totalyield + parseInt(this.summary["All reads"][readtype]["yield"]["data"]);
                readcount = readcount + parseInt(this.summary["All reads"][readtype]["read_count"]["data"]);
            }
            //console.log("yieldhistory length");
            //console.log("test" + this.livedata.live_read_count);
            if (this.livedata.yield_history.length > 1) {
                this.livedata.scalingfactor = (totalyield / readcount) / (this.livedata.yield_history[this.livedata.yield_history.length - 1][1] / this.livedata.live_read_count);
            }
        }
    }
};
