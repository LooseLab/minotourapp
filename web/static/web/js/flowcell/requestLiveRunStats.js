function parseporehist(descriptions, counts) {

    var results = [];
    var colorlookup = [];

    descriptions = JSON.parse(descriptions);

    for (var thing in descriptions["groups"]) {

        if (descriptions["groups"].hasOwnProperty(thing)) {

            if (descriptions["groups"][thing].hasOwnProperty("style")) {

                for (var state in descriptions["groups"][thing]["states"]) {

                    colorlookup[descriptions["groups"][thing]["states"][state]["name"]] = descriptions["groups"][thing]["states"][state]["style"]["colour"];
                }
            }
        }
    }

    for (var pore in counts) {

        results.push({"name": pore, "color": "#" + colorlookup[pore], "data": counts[pore]});
    }

    return results
}

function updatePoreStats(poreShizzle, livedata) {

    var returndata = parseporehist(livedata.colours_string, livedata.pore_history);

    while (poreShizzle.series.length > 0)
        poreShizzle.series[0].remove(true);

    for (var i = 0; i < returndata.length; i++) {

        var seriesdata = returndata[i];

        poreShizzle.addSeries(seriesdata);
    }
}

function updateLiveCumuYield(liveCumuYield, liveOccupancy, liveInStrand, liveTemperature, liveVoltage, livedata) {

    if (liveCumuYield.series.length < 1) {

        liveCumuYield.addSeries({
            data: livedata.yield_history
        });
        liveCumuYield.series[0].update({name: "Events"}, false);
    } else {

        liveCumuYield.series[0].setData(livedata.yield_history);
        liveCumuYield.series[0].update({name: "Events"}, false);
    }

    liveCumuYield.redraw();
    liveCumuYield.reflow();

    if (liveOccupancy.series.length < 1) {
        liveOccupancy.addSeries({
            data: livedata.percentage
        });
        liveOccupancy.series[0].update({name: "% Occupancy"}, false);
    } else {
        liveOccupancy.series[0].setData(livedata.percentage);
        liveOccupancy.series[0].update({name: "% Occupancy"}, false);
    }

    liveOccupancy.redraw();
    liveOccupancy.reflow();

    if (liveInStrand.series.length < 2) {
        liveInStrand.addSeries({data: livedata.strand});
        liveInStrand.series[0].update({name: "In Strand"}, false);
        liveInStrand.addSeries({data: livedata.good_single});
        liveInStrand.series[1].update({name: "Single Pore"}, false);
        liveInStrand.addSeries({data: livedata.adapter});
        liveInStrand.series[2].update({name: "Adapter"}, false);
        liveInStrand.addSeries({data: livedata.pore});
        liveInStrand.series[3].update({name: "Pore"}, false);

    } else {
        liveInStrand.series[0].setData(livedata.strand);
        liveInStrand.series[0].update({name: "In Strand"}, false);
        liveInStrand.series[1].setData(livedata.good_single);
        liveInStrand.series[1].update({name: "Single Pore"}, false);
        liveInStrand.series[2].setData(livedata.adapter);
        liveInStrand.series[2].update({name: "Adapter"}, false);
        liveInStrand.series[3].setData(livedata.pore);
        liveInStrand.series[3].update({name: "Pore"}, false);
    }

    liveInStrand.redraw();
    liveInStrand.reflow();

    if (liveTemperature.series.length < 1) {
        liveTemperature.addSeries({data: livedata.asictemp});
        liveTemperature.addSeries({data: livedata.heatsinktemp});
    } else {
        liveTemperature.series[0].setData(livedata.asictemp);
        liveTemperature.series[1].setData(livedata.heatsinktemp);
        liveTemperature.series[1].setData(livedata.heatsinktemp);
        liveTemperature.series[1].setData(livedata.heatsinktemp);
        liveTemperature.series[1].setData(livedata.heatsinktemp);
    }

    liveTemperature.series[0].update({name: "Asic Temp"}, false);
    liveTemperature.series[1].update({name: "HeatSink Temp"}, false);

    liveTemperature.redraw();
    liveTemperature.reflow();

    if (liveVoltage.series.length < 1) {
        liveVoltage.addSeries({data: livedata.voltage});
    } else {
        liveVoltage.series[0].setData(livedata.voltage);
    }

    liveVoltage.series[0].update({name: "Voltage"}, false);

    liveVoltage.redraw();
    liveVoltage.reflow();

}

function tohistogram(readeventcountweightedhist, readeventcountweightedhistbinwidth, totalyield) {

    var results = [];
    var categories = [];

    readeventcountweightedhist = readeventcountweightedhist.replace(/u(?=[^:]+')/g, "").replace(/'/g, "");
    readeventcountweightedhist = JSON.parse(readeventcountweightedhist);

    var totalyield = 0;
    for (var i = 0; i < readeventcountweightedhist.length; i++) {
        totalyield += readeventcountweightedhist[i];
    }

    var n50count = 0;
    var n50index = 0;
    var check = 0;

    for (i in readeventcountweightedhist) {

        n50count += parseInt(readeventcountweightedhist[i]);
        if (n50count >= (parseInt(totalyield) / 2)) {
            //console.log("n50",(i+1)*readeventcountweightedhistbinwidth, n50count);
            check += 1;
        }

        var category = String((parseInt(i)) * readeventcountweightedhistbinwidth) + " - " + String((parseInt(i) + 1) * readeventcountweightedhistbinwidth) + " ev";

        categories.push(category);

        if (check == 1) {
            n50index = i;
            results.push({"name": category, "y": parseInt(readeventcountweightedhist[i]), "color": "red"});
            check += 1;
        } else {
            results.push({"name": category, "y": parseInt(readeventcountweightedhist[i]), "color": "blue"});
        }
    }

    categories.push(">> max ev");

    var missed = 0;

    results.push({"name": ">> max ev", "y": missed});

    return [results, categories, n50index];
}

function updateLiveHistogram(liveHistogram, data) {

    var returndata = tohistogram(data[data.length - 1].minKNOW_histogram_values, data[data.length - 1].minKNOW_histogram_bin_width, data[data.length - 1].event_yield);

    var N50 = parseInt(returndata[2]);

    liveHistogram.series[0].setData(returndata[0]);

    liveHistogram.xAxis[0].setCategories(returndata[1]);

    liveHistogram.xAxis[0].removePlotBand("plot-band-1");

    liveHistogram.xAxis[0].addPlotBand({
        from: N50 - 0.5,
        to: N50 + 0.5,
        color: "#FCFFC5",
        id: "plot-band-1",
    });

    liveHistogram.xAxis[0].removePlotBand("plot-band-2");

    liveHistogram.xAxis[0].addPlotBand({
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

    liveHistogram.xAxis[0].removePlotBand("plot-band-3");

    liveHistogram.xAxis[0].addPlotBand({
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

    liveHistogram.reflow();
}

function calculatereadtoeventscaling() { // TODO is it still being used?

    var totalyield = 0;
    var readcount = 0;

    if (this.summary !== null) {

        for (var readtype in this.summary["All reads"]) {
            //console.log(this.summary["All reads"][readtype]);
            totalyield = totalyield + parseInt(this.summary["All reads"][readtype]["yield"]["data"]);
            readcount = readcount + parseInt(this.summary["All reads"][readtype]["read_count"]["data"]);
        }

        if (this.livedata.yield_history.length > 1) {
            this.livedata.scalingfactor = (totalyield / readcount) / (this.livedata.yield_history[this.livedata.yield_history.length - 1][1] / this.livedata.live_read_count);
        }
    }
}

function requestLiveRunStats(id) {

    console.log('This is lastread from flowcellcontroller: ' + flowcellController.lastread);

    // var url_livestats = '/api/v1/flowcells/' + id + '/runstats/' + flowcell_controller.lastread;
    var url_livestats = '/api/v1/flowcells/' + id + '/runstats/0';

    $.get(url_livestats, (function (result) {

        console.log('>>>> result');

        console.log(result);

        console.log('<<<< result');

        var liveHistogram = makeLiveHistogram(
            "live-histogram",
            "Histogram of Read Lengths (Events)",
            "test"
        );

        var liveCumuYield = makeLiveChart(
            "chart-yield",
            "Yield Over Time",
            "cumulative yield"
        );

        var liveInStrand = makeLiveChart(
            "live-strand",
            "in strand counts",
            "number of pores in strand/single"
        );

        var liveOccupancy = makeLiveChart(
            "live-occupancy",
            "% Occupancy Over Time",
            "% Occupancy"
        );

        var liveTemperature = makeLiveChart(
            "live-temperature",
            "Temperature Over Time",
            "°Celcius"
        );

        var liveVoltage = makeLiveChart(
            "live-voltage",
            "Voltage Over Time",
            "mV"
        );

        var poreShizzle = makeAreaPlot(
            "poreshizzle",
            "Pore States".toUpperCase(),
            "Pore States".toUpperCase()
        );

        var data = result["data"];

        var live_data = new Array();

        live_data.colours_string = result['minKNOW_colours_string'];
        live_data.voltage = new Array();
        live_data.asictemp = new Array();
        live_data.heatsinktemp = new Array();
        live_data.strand = new Array();
        live_data.adapter = new Array();
        live_data.good_single = new Array();
        live_data.pore = new Array();
        live_data.currpercentage = null;
        live_data.currstrand = null;
        live_data.percentage = new Array();
        live_data.yield_history = new Array();
        live_data.meanratio_history = new Array();
        live_data.instrand_history = new Array();
        live_data.openpore_history = new Array();
        live_data.pore_history = new Array();

        var myStringArray = ["above", "adapter", "below", "good_single", "strand", "inrange", "multiple", "pending_mux_change", "saturated", "unavailable", "unblocking", "unclassified", "unknown", "zero", "pore", "no_pore"];

        var arrayLength = myStringArray.length;

        for (var i = 0; i < arrayLength; i++) {
            live_data.pore_history[myStringArray[i]] = new Array();
        }

        if (data.length > 0) {

            data.sort(function (a, b) {
                var a_date = new Date(a.sample_time);
                var b_date = new Date(b.sample_time);

                return a_date - b_date;
            });

            for (var i = 0; i < data.length; i++) {

                timestamp = new Date(data[i].sample_time).getTime();

                live_data.live_read_count = data[i].minKNOW_read_count;
                live_data.voltage.push([timestamp, data[i].voltage_value]);
                live_data.asictemp.push([timestamp, data[i].asic_temp]);
                live_data.heatsinktemp.push([timestamp, data[i].heat_sink_temp]);
                live_data.strand.push([timestamp, data[i].strand]);
                live_data.adapter.push([timestamp, data[i].adapter]);
                live_data.pore.push([timestamp, data[i].pore]);
                live_data.good_single.push([timestamp, data[i].good_single]);
                live_data.currpercentage = data[i].occupancy;
                live_data.currstrand = data[i].strand;
                live_data.percentage.push([timestamp, data[i].occupancy]);
                live_data.yield_history.push([timestamp, data[i].event_yield]);
                live_data.meanratio_history.push([timestamp, data[i].mean_ratio]);
                live_data.instrand_history.push([timestamp, data[i].in_strand]);
                live_data.openpore_history.push([timestamp, parseInt(data[i].open_pore)]);


                for (var j = 0; j < arrayLength; j++) {

                    if (isNaN(data[i][myStringArray[j]])) {
                        live_data.pore_history[myStringArray[j]].push([timestamp, 0]);
                    } else {
                        live_data.pore_history[myStringArray[j]].push([timestamp, parseInt(data[i][myStringArray[j]])]);
                    }
                }
            }

            // calculatereadtoeventscaling();

            updateLiveHistogram(liveHistogram, data);
            updateLiveCumuYield(liveCumuYield, liveOccupancy, liveInStrand, liveTemperature, liveVoltage, live_data);
            updatePoreStats(poreShizzle, live_data);

        }

    }).bind(this));

};
