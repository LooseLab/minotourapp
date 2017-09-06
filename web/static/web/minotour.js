Highcharts.setOptions({
    plotOptions: {
        series: {
            animation: false,
            turboThreshold: 0
        }
    }
});

var NUMBER_SECONDS_IN_A_MINUTE = 60;

var MINOTOUR_VERSION = 0.5;

function check_minotour_version() {
    $.getJSON("http://www.nottingham.ac.uk/~plzloose/minoTourhome/message.php?callback=?", function (result) {

        $.each(result, function (key, value) {
            //checking version info.
            if (key == "version") {
                if (value == MINOTOUR_VERSION) {
                    $("#newstarget").html("You are running the most recent version of minoTour - version " + value + ".<br>");
                }
                else if (value < MINOTOUR_VERSION) {
                    $("#newstarget").html("You appear to be in the fortunate position of running a future version of the minoTour web application " + value + ". If you have modified the code yourself - great. If not then there might be an issue somewhere!.<br>");
                }
                else if (value > MINOTOUR_VERSION) {
                    $("#newstarget").html("You are running an outdated version of the minoTour web application. The most recent version of minoTour is version " + value + ".<br>" + "Instructions for upgrading will be posted below.<br>");
                }

            } else if (key.substring(0, 7) == "message") {
                $("#newstarget").append(value + "<br>");
            }
        });
    });
}

function check_user_runs() {
    var url = "/api/v1/runs/";

    $.getJSON(url, function (data) {
        var items = [];
        $.each(data, function (key, val) {
            items.push("<li id=" + key + "><a href='/web/private/runs/" + val.id + "'>" + val.run_name + "</a></li>");
        });

        var text = "<span>You have " + data.length + " minION runs available to view.</span>";
        text += "<ul>";
        text += items.join("");
        text += "</ul>";

        $("#div-user-runs").html(text);
    });
}


function makeChart(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: "column",
            animation: Highcharts.svg, // don"t animate in old IE
            marginRight: 10
        },
        title: {
            text: chartTitle
        },
        yAxis: {
            title: {
                text: yAxisTitle
            },
            plotLines: [{
                value: 0,
                width: 1,
                color: "#808080"
            }]
        },
        legend: {
            enabled: true
        },
        exporting: {
            enabled: false
        },
        series: []
    });

    return chart;
}


function makeChart2(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: "spline",
            marginRight: 10,
            animation: false,
            zoomType: "x"
        },
        title: {
            text: chartTitle
        },
        xAxis: {
            type: "datetime",
            tickPixelInterval: 150
        },
        yAxis: {
            title: {
                text: yAxisTitle
            },
            plotLines: [{
                value: 0,
                width: 1,
                color: "#808080"
            }]
        },
        legend: {
            enabled: true
        },
        exporting: {
            enabled: false
        }
    });

    return chart;
}


function makeLiveHistogram(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: 'column',
            //marginRight: 10,
            animation: false,
            zoomType: 'x'
        },
        title: {
            text: chartTitle,
        },
        xAxis: {
            categories: [],
        },
        yAxis: {
            title: {
                text: 'Total Event Length'
            }
        },
        credits: {
            enabled: false
        },
        series: [{
            name: 'Read Histogram',
            //data: this.datain
        }]
    });
    return chart;

}

function makeAreaPlot(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.stockChart(divName, {
        chart: {
            renderTo: 'container-porehist' + this.title,
            type: 'area',
            //type: 'spline',
            height: 350,
        },
        title: {
            text: chartTitle,
        },
        xAxis: {
            range: 1 * 360 * 1000, //set range to last hour of data
        },
        rangeSelector: {
            enabled: false
        },
        yAxis: {
            //max: 512,
            endOnTick: false,
            title: {
                text: 'Channel Classifications'
            }
        },
        legend: {
            enabled: true
        },
        plotOptions: {
            area: {
                stacking: 'percent',
            },
            series: {
                showInNavigator: true,
                dataLabels: {
                    enabled: false,
                    formatter: function () {
                        return this.y;
                    }
                }
            }
        },
        credits: {
            enabled: false
        },
        series: [],
    });
    return chart;
}

function makeYieldProjection(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: 'spline',
            zoomType: 'x',
            height: 350,
        },
        title: {
            text: chartTitle,
        },
        xAxis: {
            type: 'datetime',
            tickPixelInterval: 150
        },
        yAxis: {
            title: {
                text: 'Cumulative Predicted Bases'
            },
            plotLines: [{
                value: 0,
                width: 1,
                color: '#808080'
            }],
            min: 0,
        },
        credits: {
            enabled: false
        },
        series: [{
            name: 'Original Data',
            data: []
        },
            {
                name: 'Current Projected Data',
                dashStyle: 'longdash',
                data: []
            },
            {
                name: 'First Hour Projected Data',
                dashStyle: 'shortdash',
                data: []
            },
            {
                name: 'Ideal Results',
                dashStyle: 'Dot',
                data: []
            }
        ]

    });
    return chart;
}

function makeLiveChart(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.stockChart(divName, {
        chart: {
            type: 'spline',
            zoomType: 'xy'
        }
        ,
        rangeSelector: {
            enabled: true,
            buttons:
                [{
                    type: 'minute',
                    count: 1,
                    text: '1min'
                }, {
                    type: 'minute',
                    count: 5,
                    text: '5min'
                }, {
                    type: 'minute',
                    count: 30,
                    text: '1/2hr'
                }, {
                    type: 'minute',
                    count: 60,
                    text: '1hr'
                }, {
                    type: 'day',
                    count: 0.5,
                    text: '12hrs'
                }, {
                    type: 'day',
                    count: 1,
                    text: '1day'
                }, {
                    type: 'all',
                    text: 'All'
                }]
        }
        ,
        title: {
            text: chartTitle,
        }
        ,
        xAxis: {
            type: 'datetime',
            tickPixelInterval: 150
        }
        ,
        yAxis: {
            title: {
                text: yAxisTitle
            }
            ,
            plotLines: [{
                value: 0,
                width: 1,
                color: '#808080'
            }],
            min: 0,
        }
        ,
        credits: {
            enabled: false
        }
        ,
        series: []
    });
    return chart;
};

function MonitorAPP() {
    var livedata = new Array();
    var self = this;

    this.init = function () {
        this.requestData();

        setInterval(function () {
            this.requestData();
        }.bind(this), 10000);
    };

    this.updatecounters = function (live) {
        var thing = document.getElementById("livenum");
        thing.innerHTML = live;
        var thing2 = document.getElementById("livenumruns");
        thing2.innerHTML = "You have " + live + " live runs.";
    };

    this.requestData = function () {
        var url_run = "/api/v1/currentruns/";

        $.get(url_run, function (data) {
            //console.log(data);
            self.livedata = data;
            self.updatecounters(self.livedata.length);
            //self.barcodes = data.barcodes.sort();
            //self.updateBarcodeNavTab();
        }.bind(this));

        //self.updatecounters(self.livedata.length + "/" + "0");
        //self.requestSummaryByMinuteData(self.id);
        //self.requestSummaryData(self.id);
    };
}

function MinotourApp() {
    this.chart_reads_called = null;
    this.chart_yield = null;
    this.chart_average_read_length = null;
    this.chart_maximum_read_length = null;
    this.average_read_lengths_overtime = null;
    this.chart_cumulative_number_reads_overtime = null;
    this.chartSequencingRate = null;
    this.chartHistogramReadLength = null;
    this.chartHistogramBasesSequencedByReadLength = null;

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
    }
    ;
    this.livedata.minIONname = null;
    this.livedata.colours_string = null;
    this.livedata.scalingfactor = 0;

    this.makeChart = makeChart;
    this.makeChart2 = makeChart2;
    this.makeLiveHistogram = makeLiveHistogram;
    this.makeYieldProjection = makeYieldProjection;
    this.makeLiveChart = makeLiveChart;
    this.makeAreaPlot = makeAreaPlot;
    this.lastread = 0;
    this.needtoupdatecharts = false;

    var self = this;

    this.init = function () {
        /*
        * Add event listeners
        */

        document.getElementById("panel-live-data").style.display = "none";
        document.getElementById("panel-run-summary-data").style.display = "none";
        document.getElementById("panel-basecalled-data").style.display = "block";

        document.getElementById("run-summary-data").onclick = function (e) {
            document.getElementById("panel-basecalled-data").style.display = "none";
            document.getElementById("nav-basecalled-data").parentNode.classList.remove("active");
            document.getElementById("panel-live-data").style.display = "none";
            document.getElementById("nav-live-data").parentNode.classList.remove("active");
            document.getElementById("panel-run-summary-data").style.display = "block";
            document.getElementById("run-summary-data").parentNode.classList.add("active");
        };

        document.getElementById("nav-live-data").onclick = function (e) {
            document.getElementById("panel-basecalled-data").style.display = "none";
            document.getElementById("nav-basecalled-data").parentNode.classList.remove("active");
            document.getElementById("panel-live-data").style.display = "block";
            document.getElementById("nav-live-data").parentNode.classList.add("active");
            document.getElementById("panel-run-summary-data").style.display = "none";
            document.getElementById("run-summary-data").parentNode.classList.remove("active");
        };

        document.getElementById("nav-basecalled-data").onclick = function (e) {
            document.getElementById("panel-basecalled-data").style.display = "block";
            document.getElementById("nav-basecalled-data").parentNode.classList.add("active");
            document.getElementById("panel-live-data").style.display = "none";
            document.getElementById("nav-live-data").parentNode.classList.remove("active");
            document.getElementById("panel-run-summary-data").style.display = "none";
            document.getElementById("run-summary-data").parentNode.classList.remove("active");
        };

        this.PoreShizzle = this.makeAreaPlot(
            'poreshizzle',
            'Pore States'.toUpperCase(),
            'Pore States'.toUpperCase()
        );

        this.chart_reads_called = this.makeChart(
            "reads-called",
            "reads called".toUpperCase(),
            "number of reads called".toUpperCase()
        );

        this.chart_yield = this.makeChart(
            "yield",
            "yield".toUpperCase(),
            "yield".toUpperCase()
        );

        this.chart_average_read_length = this.makeChart(
            "average-read-length",
            "average read length".toUpperCase(),
            "average read length".toUpperCase()
        );

        this.chart_maximum_read_length = this.makeChart(
            "maximum-read-length",
            "maximum read length".toUpperCase(),
            "maximum read length".toUpperCase()
        );

        this.average_read_lengths_overtime = this.makeChart2(
            "average-read-lengths-overtime",
            "average read length over time".toUpperCase(),
            "average read length".toUpperCase()
        );

        this.chart_cumulative_number_reads_overtime = this.makeChart2(
            "cumulative-number-reads-overtime",
            "cumulative reads".toUpperCase(),
            "cumulative reads".toUpperCase()
        );

        this.chartSequencingRate = this.makeChart2(
            "sequencing-rate",
            "sequencing rate".toUpperCase(),
            "bases/second".toUpperCase()
        );

        this.chartSequencingSpeed = this.makeChart2(
            "sequencing-speed",
            "sequencing speed".toUpperCase(),
            "bases/channel/second".toUpperCase()
        );

        this.chartHistogramReadLength = this.makeChart2(
            "histogram-read-lengths",
            "Histogram of Read Lengths".toUpperCase(),
            "Number of reads".toUpperCase()
        );

        this.chartHistogramBasesSequencedByReadLength = this.makeChart2(
            "histogram-bases-sequenced-by-read-length",
            "Histogram of Bases Sequenced by Read Length".toUpperCase(),
            "Number of bases".toUpperCase()
        );

        this.LiveHistogram = this.makeLiveHistogram(
            'live-histogram',
            'Histogram of Read Lengths (Events)',
            'test'
        );

        this.LiveYield = this.makeYieldProjection(
            'yield-projection',
            'Yield Projection',
            'Yield Projection'
        );

        this.LiveCumuYield = this.makeLiveChart(
            'chart-yield',
            'Yield Over Time',
            'cumulative yield'
        );

        this.LiveInStrand = this.makeLiveChart(
            'live-strand',
            'in strand counts',
            'number of pores in strand/single'
        );

        this.LiveOccupancy = this.makeLiveChart(
            'live-occupancy',
            '% Occupancy Over Time',
            '% Occupancy'
        );

        this.LiveTemperature = this.makeLiveChart(
            'live-temperature',
            'Temperature Over Time',
            '°Celcius'
        );
        this.LiveVoltage = this.makeLiveChart(
            'live-voltage',
            'Voltage Over Time',
            'mV'
        );
        this.LivePoreState = this.makeLiveChart(
            'live-porestate',
            'Pore State Currents',
            'Current pA'
        );
        this.LiveCurrentRatio = this.makeLiveChart(
            'live-currentratio',
            'Current Ratio In Strand/Open Pore',
            'Current Ratio'
        );

        self.id = document.getElementById("run-id").innerText;
        self.selectedBarcode = "All reads";

        this.requestData();

        setInterval(function () {
            this.requestData();
        }.bind(this), 10000);
    };

    this.numberWithCommas = function (x) {
        return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
    }

    this.updatetext = function (livedata) {
        document.getElementById('MinIONName').innerHTML = self.livedata.minIONname;
        document.getElementById('asicid').innerHTML = self.livedata.asicid;
        document.getElementById('MinKNOWVersion').innerHTML = livedata.minKNOW_version;
        document.getElementById('RunID').innerHTML = livedata.run_id;
        document.getElementById('RunName').innerHTML = livedata.run_name;
        document.getElementById('SampleName').innerHTML = livedata.sample_name;
        //document.getElementById('ComputerName').innerHTML = livedata.minKNOW_computer;
        var starttime = new Date(livedata.start_time);
        document.getElementById('RunStartTime').innerHTML = starttime;
        document.getElementById('FlowCellID').innerHTML = livedata.minKNOW_flow_cell_id;
        if (livedata.active == true) {
            document.getElementById('CurrentlySequencing').innerHTML = '<i class="fa fa-check" aria-hidden="true"></i>';
        } else {
            document.getElementById('CurrentlySequencing').innerHTML = '<i class="fa fa-times" aria-hidden="true"></i>';
        }
        if (livedata.barcodes.length > 2) {
            document.getElementById('Barcoded').innerHTML = '<i class="fa fa-check" aria-hidden="true"></i>';
        } else {
            document.getElementById('Barcoded').innerHTML = '<i class="fa fa-times" aria-hidden="true"></i>';
        }


    }

    /*
     * Each click on a barcode tab fires this function
     * that calls requestData and update all charts
     */
    this.updateChartsBasedOnBarcode = function (event) {
        self.selectedBarcode = event.target.innerText;
        self.requestData();
    };

    /*
     * Updates the list of barcodes tab and attach
     * click event to function updateChartsBasedOnBarcode
     */
    this.updateBarcodeNavTab = function () {
        var ul = document.getElementById("nav-tabs-barcodes");

        ul.innerHTML = "";

        var sortedBarcodes = this.barcodes;

        for (var i = 0; i < sortedBarcodes.length; i++) {
            var li = document.createElement("li");
            var a = document.createElement("a");
            a.onclick = this.updateChartsBasedOnBarcode;
            a.href = "#";
            a.text = sortedBarcodes[i];

            if (sortedBarcodes[i] === self.selectedBarcode) {
                li.classList.add("active");
            }

            li.appendChild(a);
            ul.appendChild(li);
        }
    };

    this.updateReadsColumnBasedChart = function (chart, field) {
        var summaries = this.summary;

        var series = [];

        // Always include all reads
        data = [];

        for (var readtype of Object.keys(summaries["All reads"])) {
            data.push(summaries["All reads"][readtype][field]["data"][0]);
        }

        serie = {
            "name": "All reads",
            "data": data
        };

        series.push(serie);

        // Include specific barcode if selected
        if (self.selectedBarcode !== "All reads") {

            data = [];

            for (var readtype of Object.keys(summaries[self.selectedBarcode])) {
                data.push(summaries[self.selectedBarcode][readtype][field]["data"][0]);
            }

            serie = {
                "name": self.selectedBarcode,
                "data": data
            };

            series.push(serie);

        }

        chart.colorCounter = 2;
        chart.symbolCounter = 0;

        var chartSeriesLength = (chart.series ? chart.series.length : 0);

        for (var i = 0; i < series.length; i++) {
            if (i <= (chartSeriesLength - 1)) {
                chart.series[i].setData(series[i].data);
                chart.series[i].update({
                    name: series[i].name
                });
            } else {
                chart.addSeries(series[i]);
            }
        }

        chartSeriesLength = (chart.series ? chart.series.length : 0);

        while (chartSeriesLength > series.length) {
            chart.series[(chartSeriesLength - 1)].remove();
            chartSeriesLength = (chart.series ? chart.series.length : 0);
        }
    };

    this.updateAverageReadLengthOverTimeChart = function () {

        var chart = self.average_read_lengths_overtime;
        var selectedBarcode = self.selectedBarcode;

        while (chart.series.length > 0) {
            chart.series[0].remove();
        }

        if (selectedBarcode !== "All reads") {
            var summaries = {};
            summaries["All reads"] = {};
            summaries[selectedBarcode] = {};

        } else {
            var summaries = {
                "All reads": {}
            };

        }

        for (var i = 0; i < self.summaryByMinute.length; i++) {

            var average_read_length = self.summaryByMinute[i].total_length / self.summaryByMinute[i].read_count;
            var sample_time = new Date(self.summaryByMinute[i].sample_time);

            var point = {
                x: sample_time,
                y: average_read_length
            }

            if (self.summaryByMinute[i].barcodename === "All reads") {
                if (summaries["All reads"][self.summaryByMinute[i].typename] === undefined) {
                    summaries["All reads"][self.summaryByMinute[i].typename] = [];
                }

                summaries["All reads"][self.summaryByMinute[i].typename].push(point);
            }

            if (self.summaryByMinute[i].barcodename === selectedBarcode && selectedBarcode !== "All reads") {
                if (summaries[selectedBarcode][self.summaryByMinute[i].typename] === undefined) {
                    summaries[selectedBarcode][self.summaryByMinute[i].typename] = [];
                }

                summaries[selectedBarcode][self.summaryByMinute[i].typename].push(point);
            }

        }

        for (var barcode in summaries) {
            for (var readtype in summaries[barcode]) {
                chart.addSeries({name: barcode + " - " + readtype, data: summaries[barcode][readtype]});
            }
        }
    };

    this.updateCumulativeNumberOfReadsOverTimeChart = function () {
        var chart = self.chart_cumulative_number_reads_overtime;
        var selectedBarcode = self.selectedBarcode;

        while (chart.series.length > 0) {
            chart.series[0].remove();
        }

        if (selectedBarcode !== "All reads") {
            var summaries = {};
            summaries["All reads"] = {};
            summaries[selectedBarcode] = {};

        } else {
            var summaries = {
                "All reads": {}
            };

        }

        for (var i = 0; i < self.summaryByMinute.length; i++) {

            if (self.summaryByMinute[i].barcodename === "All reads") {
                if (summaries["All reads"][self.summaryByMinute[i].typename] === undefined) {
                    summaries["All reads"][self.summaryByMinute[i].typename] = {
                        "lastCumulativeReadCount": 0,
                        "data": []
                    };
                }

                var cumulativeReadCount = self.summaryByMinute[i].read_count + summaries["All reads"][self.summaryByMinute[i].typename].lastCumulativeReadCount;
                var sample_time = new Date(self.summaryByMinute[i].sample_time);

                var point = {
                    x: sample_time,
                    y: cumulativeReadCount
                }

                summaries["All reads"][self.summaryByMinute[i].typename].lastCumulativeReadCount = cumulativeReadCount;
                summaries["All reads"][self.summaryByMinute[i].typename].data.push(point);
            }

            if (self.summaryByMinute[i].barcodename === selectedBarcode && selectedBarcode !== "All reads") {
                if (summaries[selectedBarcode][self.summaryByMinute[i].typename] === undefined) {
                    summaries[selectedBarcode][self.summaryByMinute[i].typename] = {
                        "lastCumulativeReadCount": 0,
                        "data": []
                    };
                }

                var cumulativeReadCount = self.summaryByMinute[i].read_count + summaries[selectedBarcode][self.summaryByMinute[i].typename].lastCumulativeReadCount;
                var sample_time = new Date(self.summaryByMinute[i].sample_time);

                var point = {
                    x: sample_time,
                    y: cumulativeReadCount
                }

                summaries[selectedBarcode][self.summaryByMinute[i].typename].lastCumulativeReadCount = cumulativeReadCount;
                summaries[selectedBarcode][self.summaryByMinute[i].typename].data.push(point);
            }

        }

        for (var barcode in summaries) {
            for (var readtype in summaries[barcode]) {
                chart.addSeries({name: barcode + " - " + readtype, data: summaries[barcode][readtype]["data"]});
            }
        }
    };

    this.updateSequencingRateChart = function () {
        var chart = self.chartSequencingRate;
        var selectedBarcode = self.selectedBarcode;

        // Remove previous series
        while (chart.series.length > 0) {
            chart.series[0].remove();
        }

        for (var barcode of Object.keys(self.summaryByMinute2)) {

            if (barcode === 'All reads' || barcode === self.selectedBarcode) {
                for (var typeName of Object.keys(self.summaryByMinute2[barcode])) {

                    chart.addSeries({
                        name: barcode + " - " + typeName,
                        data: self.summaryByMinute2[barcode][typeName]["sequencingRate"]
                    });

                }
            }
        }

    };

    this.updateSequencingSpeedChart = function () {
        var chart = self.chartSequencingSpeed;
        var selectedBarcode = self.selectedBarcode;

        // Remove previous series
        while (chart.series.length > 0) {
            chart.series[0].remove();
        }

        for (var barcode of Object.keys(self.summaryByMinute2)) {

            if (barcode === 'All reads' || barcode === self.selectedBarcode) {
                for (var typeName of Object.keys(self.summaryByMinute2[barcode])) {

                    chart.addSeries({
                        name: barcode + " - " + typeName,
                        data: self.summaryByMinute2[barcode][typeName]["sequencingSpeed"]
                    });

                }
            }
        }

    };

    this.updateHistogramReadLengthChart = function () {

        var chart = self.chartHistogramReadLength;

        // Remove previous series
        while (chart.series.length > 0) {
            chart.series[0].remove();
        }

        //console.log('Inside function updateHistogramReadLengthChart.');
        //console.log(Object.keys(self.histogramSummary));

        for (var barcode_name of Object.keys(self.histogramSummary)) {

            //console.log('barcode_name: '+ barcode_name + ', self.selectedBarcode: ' + self.selectedBarcode);
            //console.log(barcode_name === self.selectedBarcode);

            //if (barcode_name === 'All reads' || barcode_name === self.selectedBarcode) {
            if (barcode_name === self.selectedBarcode) {

                for (var typeName of Object.keys(self.histogramSummary[barcode_name])) {

                    //console.log('typeName: ' + typeName);
                    //console.log(self.histogramSummary[barcode_name][typeName]['bin_width']);
                    //console.log(self.histogramSummary[barcode_name][typeName]["read_count"]);

                    chart.update({
                        chart: {
                            type: 'column'
                        },
                        xAxis: {
                            type: 'category',
                            categories: self.histogramSummary[barcode_name][typeName]['bin_width']
                        }
                    });

                    chart.addSeries({
                        name: barcode_name + " - " + typeName,
                        data: self.histogramSummary[barcode_name][typeName]["read_count"]
                    });

                }
            }
        }

    };

    this.updateHistogramBasesSequencedReadLengthChart = function () {

        var chart = self.chartHistogramBasesSequencedByReadLength;

        // Remove previous series
        while (chart.series.length > 0) {
            chart.series[0].remove();
        }

        //console.log('Inside function updateHistogramReadLengthChart.');
        //console.log(Object.keys(self.histogramSummary));

        for (var barcode_name of Object.keys(self.histogramSummary)) {

            //console.log('barcode_name: '+ barcode_name + ', self.selectedBarcode: ' + self.selectedBarcode);
            //console.log(barcode_name === self.selectedBarcode);

            //if (barcode_name === 'All reads' || barcode_name === self.selectedBarcode) {
            if (barcode_name === self.selectedBarcode) {

                for (var typeName of Object.keys(self.histogramSummary[barcode_name])) {

                    //console.log('typeName: ' + typeName);
                    //console.log(self.histogramSummary[barcode_name][typeName]['bin_width']);
                    //console.log(self.histogramSummary[barcode_name][typeName]["read_count"]);

                    chart.update({
                        chart: {
                            type: 'column'
                        },
                        xAxis: {
                            type: 'category',
                            categories: self.histogramSummary[barcode_name][typeName]['bin_width']
                        }
                    });

                    chart.addSeries({
                        name: barcode_name + " - " + typeName,
                        data: self.histogramSummary[barcode_name][typeName]["read_length"]
                    });

                }
            }
        }

    };

    this.updateSummaryByMinuteBasedCharts = function () {
        self.updateAverageReadLengthOverTimeChart()
        self.updateCumulativeNumberOfReadsOverTimeChart();
        self.updateSequencingRateChart();
        self.updateSequencingSpeedChart();
    };

    this.updateLiveCumuYield = function () {
        //console.log(self.LiveCumuYield);
        if (self.LiveCumuYield.series.length < 1){
            self.LiveCumuYield.addSeries({
                data: self.livedata.yield_history });
            self.LiveCumuYield.series[0].update({name: "Events"}, false);
        }else {
            self.LiveCumuYield.series[0].setData(self.livedata.yield_history);
            self.LiveCumuYield.series[0].update({name: "Events"}, false);
        }
        self.LiveCumuYield.redraw();
        self.LiveCumuYield.reflow();
        if (self.LiveOccupancy.series.length < 1){
            self.LiveOccupancy.addSeries({
                data: self.livedata.percentage });
            self.LiveOccupancy.series[0].update({name: "% Occupancy"}, false);
        }else {
            self.LiveOccupancy.series[0].setData(self.livedata.percentage);
            self.LiveOccupancy.series[0].update({name: "% Occupancy"}, false);
        }
        self.LiveOccupancy.redraw();
        self.LiveOccupancy.reflow();
        if (self.LiveInStrand.series.length < 2){
            self.LiveInStrand.addSeries({data:self.livedata.strand});
            self.LiveInStrand.series[0].update({name: "In Strand"}, false);
            self.LiveInStrand.addSeries({data:self.livedata.good_single});
            self.LiveInStrand.series[1].update({name: "Single Pore"}, false);

        }else {
            self.LiveInStrand.series[0].setData(self.livedata.strand);
            self.LiveInStrand.series[0].update({name: "In Strand"}, false);
            self.LiveInStrand.series[1].setData(self.livedata.good_single);
            self.LiveInStrand.series[1].update({name: "Single Pore"}, false);
        }
        self.LiveInStrand.redraw();
        self.LiveInStrand.reflow();
        if (self.LivePoreState.series.length < 1){
            self.LivePoreState.addSeries({data: self.livedata.instrand_history});
            self.LivePoreState.addSeries({data: self.livedata.openpore_history});
        }else{
            self.LivePoreState.series[0].setData(self.livedata.instrand_history);
            self.LivePoreState.series[1].setData(self.livedata.openpore_history);
        }
        self.LivePoreState.series[0].update({name: "In Strand"}, false);
        self.LivePoreState.series[1].update({name: "Open Pore"}, false);
        self.LivePoreState.redraw();
        self.LivePoreState.reflow();
        if (self.LiveCurrentRatio.series.length < 1){
            self.LiveCurrentRatio.addSeries({data: self.livedata.meanratio_history});
        }else{
            self.LiveCurrentRatio.series[0].setData(self.livedata.meanratio_history);
        }
        self.LiveCurrentRatio.series[0].update({name: "Current Ratio"}, false);
        self.LiveCurrentRatio.redraw();
        self.LiveCurrentRatio.reflow();
        if (self.LiveTemperature.series.length < 1){
            self.LiveTemperature.addSeries({data: self.livedata.asictemp});
            self.LiveTemperature.addSeries({data: self.livedata.heatsinktemp});
        }else{
            self.LiveTemperature.series[0].setData(self.livedata.asictemp);
            self.LiveTemperature.series[1].setData(self.livedata.heatsinktemp);
        }
        self.LiveTemperature.series[0].update({name: "Asic Temp"}, false);
        self.LiveTemperature.series[1].update({name: "HeatSink Temp"}, false);
        self.LiveTemperature.redraw();
        self.LiveTemperature.reflow();
        if (self.LiveVoltage.series.length < 1){
            self.LiveVoltage.addSeries({data: self.livedata.voltage});
        }else{
            self.LiveVoltage.series[0].setData(self.livedata.voltage);
        }
        self.LiveVoltage.series[0].update({name: "Voltage"}, false);
        self.LiveVoltage.redraw();
        self.LiveVoltage.reflow();

    };

    this.updateSummaryBasedCharts = function () {
        var charts = {
            "read_count": self.chart_reads_called,
            "yield": self.chart_yield,
            "average_read_length": self.chart_average_read_length,
            "max_length": self.chart_maximum_read_length
        };

        for (var prop in charts) {
            self.updateReadsColumnBasedChart(charts[prop], prop);
        }
    };

    this.updateHistogramBasedCharts = function () {
        self.updateHistogramReadLengthChart();
        self.updateHistogramBasesSequencedReadLengthChart();
    }

    this.requestSummaryByMinuteData = function (id) {
        /*
         * Request summary by minute data
         */
        var url = "/api/v1/runs/" + id + "/summarybarcodebyminute";

        $.get(url, function (data) {

            if (data.length > 0) {
                var orderedData = data.sort(function (a, b) {
                    return new Date(a.sample_time) - new Date(b.sample_time);
                });

                self.summaryByMinute = orderedData;

                /*********/

                var summaries = {};

                for (var barcode of self.barcodes) {
                    summaries[barcode] = {};
                }

                for (var i = 0; i < self.summaryByMinute.length; i++) {
                    var item = self.summaryByMinute[i];

                    if (summaries[item.barcodename][item.typename] === undefined) {
                        summaries[item.barcodename][item.typename] = {
                            "data": [],
                            "sequencingRate": [],
                            "sequencingSpeed": [],
                        };
                    }

                    var sampleTime = new Date(item.sample_time);

                    var singleData = {
                        sampleTime: sampleTime,
                        totalLength: item.total_length,
                        readCount: item.read_count,
                        maxLength: item.max_length,
                        minLength: item.min_length,
                    }

                    summaries[item.barcodename][item.typename]["sequencingRate"].push({
                        x: sampleTime,
                        y: item.total_length / NUMBER_SECONDS_IN_A_MINUTE
                    });

                    summaries[item.barcodename][item.typename]["sequencingSpeed"].push({
                        x: sampleTime,
                        y: item.total_length / item.number_active_channels / NUMBER_SECONDS_IN_A_MINUTE
                    });

                }

                self.summaryByMinute2 = summaries;

                /*********/

                // update chart - TODO split ajax call and update of self.summaryByMinute from the redrawing the charts
                self.updateSummaryByMinuteBasedCharts();
            }
        });
    };

    this.calculatereadtoeventscaling = function () {
        var totalyield = 0;
        var readcount = 0;
        for (var readtype in self.summary["All reads"]) {
            //console.log(self.summary['All reads'][readtype]);
            totalyield = totalyield + parseInt(self.summary['All reads'][readtype]['yield']['data']);
            readcount = readcount + parseInt(self.summary['All reads'][readtype]['read_count']['data']);
        }
        self.livedata.scalingfactor = (totalyield / readcount) / (self.livedata.yield_history[self.livedata.yield_history.length - 1][1] / self.livedata.live_read_count);
    };

    this.requestSummaryData = function (id) {
        /*
         * Request summary by barcode data
         */
        var url = "/api/v1/runs/" + id + "/summarybarcode";

        $.get(url, function (data) {

            if (data.length > 0) {
                var summaries = {};
                for (var i = 0; i < data.length; i++) {
                    var item = data[i];

                    if (summaries[item.barcodename] === undefined) {
                        summaries[item.barcodename] = {};
                    }

                    if (summaries[item.barcodename][item.typename] === undefined) {
                        summaries[item.barcodename][item.typename] = {
                            "read_count": null,
                            "yield": null,
                            "average_read_length": null,
                            "max_length": null
                        };

                        summaries[item.barcodename][item.typename]["read_count"] = {
                            "name": item.typename,
                            "data": [item.read_count],
                            "animation": false
                        };

                        summaries[item.barcodename][item.typename]["yield"] = {
                            "name": item.typename,
                            "data": [item.total_length],
                            "animation": false
                        };

                        summaries[item.barcodename][item.typename]["average_read_length"] = {
                            "name": item.typename,
                            "data": [item.total_length / item.read_count],
                            "animation": false
                        };

                        summaries[item.barcodename][item.typename]["max_length"] = {
                            "name": item.typename,
                            "data": [item.max_length],
                            "animation": false
                        };
                    }
                }

                self.summary = summaries;

                self.updateSummaryBasedCharts();
                self.calculatereadtoeventscaling();
            }
        });

    };

    this.requestHistogramData = function (id) {
        /*
         * Request histogram data
         */

        var url = "/api/v1/runs/" + id + "/histogramsummary";

        $.get(url, function (data) {

            if (data.length > 0) {

                var ordered_data = data.sort(function (a, b) {
                    return a.bin_width - b.bin_width;
                });

                var summaries = {};

                for (var i = 0; i < ordered_data.length; i++) {

                    var item = ordered_data[i];

                    if (summaries[item.barcode_name] === undefined) {
                        summaries[item.barcode_name] = {};
                    }

                    if (summaries[item.barcode_name][item.read_type_name] === undefined) {
                        summaries[item.barcode_name][item.read_type_name] = {
                            'bin_width': [],
                            'read_count': [],
                            'read_length': []
                        };
                    }

                    summaries[item.barcode_name][item.read_type_name]['read_count'].push(item.read_count);
                    summaries[item.barcode_name][item.read_type_name]['read_length'].push(item.read_length);
                    summaries[item.barcode_name][item.read_type_name]['bin_width'].push(parseInt(item.bin_width));
                }

                self.histogramSummary = summaries;

                self.updateHistogramBasedCharts();

            }
        });
    }

    this.parseporehist = function (descriptions, counts) {
        var results = [];
        var colors = [];
        var categories = [];
        var datam = [];
        var colorlookup = [];
        descriptions = JSON.parse(descriptions);
        //console.log(descriptions);
        for (var thing in descriptions) {
            if (descriptions.hasOwnProperty(thing)) {
                if (descriptions[thing].hasOwnProperty("style")) {
                    //console.log(descriptions[thing]["style"]["colour"]);
                    colorlookup[descriptions[thing]["name"]] = descriptions[thing]["style"]["colour"];
                }
            }
        }
        for (var pore in counts) {
            results.push({"name": pore, "color": "#" + colorlookup[pore], "data": counts[pore]})//,"color":"#121212"]});
        }
        return results
    }


    this.tohistogram = function (readeventcountweightedhist, readeventcountweightedhistbinwidth, totalyield) {
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
                //console.log('n50',(i+1)*readeventcountweightedhistbinwidth, n50count);
                check += 1;
            }
            //console.log(i);
            //console.log(parseInt(i)+1);
            var category = String((parseInt(i)) * readeventcountweightedhistbinwidth) + " - " + String((parseInt(i) + 1) * readeventcountweightedhistbinwidth) + " ev";
            categories.push(category);
            if (check == 1) {
                n50index = i;
                results.push({"name": category, "y": parseInt(readeventcountweightedhist[i]), "color": 'red'});
                check += 1;
            } else {
                results.push({"name": category, "y": parseInt(readeventcountweightedhist[i]), "color": 'blue'});
            }

            //}
        }
        categories.push(">> max ev");
        var missed = 0;
        //var missed = totalyield - readeventcountweightedhist.reduce(add,0);
        results.push({"name": ">> max ev", "y": missed});
        //console.log(n50index);
        return [results, categories, n50index];

    };

    this.geteighthours = function (data, runstart) {
        return (28800 * 1000 - (data[0][0] - (runstart * 1000)));
    };

    this.converttobases = function (data, seqspeed) {
        if (self.livedata.scalingfactor > 0) {
            //console.log("returning scaling factor");
            return self.livedata.scalingfactor;
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
            var scaleddata = [];
            for (var i = 0; i < data.length; i++) {
                scaleddata.push([data[i][0], data[i][1] * scaling]);
            }
            //console.log("returning estimated scaling factor");
            return scaleddata;
        }
    };

    this.projectdata = function (data) {
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
            return [1, meanholder / (holder.length - 1 )];
        } else if (diffholder / (holder.length - 1) < 0.999) {
            return ([0.999, meanholder / (holder.length - 1 )]);
        } else {
            return ([diffholder / (holder.length - 1), meanholder / (holder.length - 1 )]);
        }
    };

    this.projectresults = function (syntheticdata, scalingfactor, steps, difference, runstart) {
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
    };

    this.scaleyield = function (firstelement, data) {
        var results = [];
        for (var i = 0; i < data.length; i++) {
            //console.log(data[i]);
            //console.log((data[i][0]-firstelement[0])/1000);
            //console.log((data[i][1]/1000000));
            results.push((((data[i][0] - firstelement[0]) / 1000), Math.ceil((data[i][1] / 1000000))));
        }
        return results;
    };


    this.updateLiveYieldProjection = function () {
        var seqspeed = "450 b/s";
        //commented out below as it wasn't being used - not sure why - need to check
        //var timeleft = self.geteighthours(self.livedatayield.slice(-1), self.rundata.start_time);
        var firsthour = 120;
        [scalingfactor, difference] = self.projectdata(self.livedata.yield_history);
        [scalingfactor2, difference2] = self.projectdata(self.livedata.yield_history.slice(0, firsthour));
        //synthericdata needs renaming!!!
        var syntheticdata = self.livedata.yield_history;
        //console.log(scalingfactor2);
        //console.log(syntheticdata);
        newarray = self.projectresults(syntheticdata, scalingfactor, 4000, difference, new Date(self.rundata.start_time).getTime());
        //console.log(newarray);
        newarray1 = self.projectresults(syntheticdata.slice(0, firsthour), scalingfactor2, 4000, difference2, new Date(self.rundata.start_time).getTime());
        newarray2 = self.projectresults(syntheticdata.slice(0, firsthour), 1, 4000, difference2, new Date(self.rundata.start_time).getTime());
        self.LiveYield.series[0].setData(self.converttobases(self.livedata.yield_history, seqspeed));
        self.LiveYield.series[1].setData(self.converttobases(newarray, seqspeed));
        self.LiveYield.series[2].setData(self.converttobases(newarray1, seqspeed));
        self.LiveYield.series[3].setData(self.converttobases(newarray2, seqspeed));
        self.LiveYield.redraw();
        self.LiveYield.reflow()
    };

    this.updatePoreStats = function () {
        var returndata = self.parseporehist(this.livedata.colours_string, this.livedata.pore_history);
        //console.log(returndata);
        while (self.PoreShizzle.series.length > 0)
            self.PoreShizzle.series[0].remove(true);
        //self.PoreShizzle.addSeries(returndata[4]);
        //console.log(returndata[4]);
        for (var i = 0; i < returndata.length; i++) {
            //console.log(returndata[i]);
            self.PoreShizzle.addSeries(returndata[i]);
        }
    }
    this.updateLiveHistogram = function (data) {
        returndata = self.tohistogram(data[data.length - 1].minKNOW_histogram_values, data[data.length - 1].minKNOW_histogram_bin_width, data[data.length - 1].event_yield);
        self.LiveHistogram.series[0].setData(returndata[0]);
        self.LiveHistogram.xAxis[0].setCategories(returndata[1]);
        var N50 = parseInt(returndata[2]);
        self.LiveHistogram.xAxis[0].removePlotBand('plot-band-1');
        self.LiveHistogram.xAxis[0].addPlotBand({
            from: N50 - 0.5,
            to: N50 + 0.5,
            color: '#FCFFC5',
            id: 'plot-band-1',
        });
        self.LiveHistogram.xAxis[0].removePlotBand('plot-band-2');
        self.LiveHistogram.xAxis[0].addPlotBand({
            color: 'black',
            width: 2,
            dashStyle: 'longdashdot',
            value: returndata[2],
            label: {
                text: 'Estimated Read N50',
                align: 'left',
                rotation: 0,
                x: +10 // Amount of pixels the label will be repositioned according to the alignment.
            },
            id: 'plot-band-2',
        });
        self.LiveHistogram.xAxis[0].removePlotBand('plot-band-3');
        self.LiveHistogram.xAxis[0].addPlotBand({
            color: 'black',
            width: 2,
            dashStyle: 'longdashdot',
            value: (Math.floor(this.totalyield / this.readcount / this.datain2)),
            label: {
                text: 'Estimated Read Average - ' + Math.round(this.totalyield / this.readcount / 1000 * 100) / 100 + ' K events',
                align: 'left',
                rotation: 0,
                x: +10,
                y: +30, // Amount of pixels the label will be repositioned according to the alignment.
            },
            id: 'plot-band-3',
        });
        self.LiveHistogram.reflow();
    };

    this.requestRunDetails = function (id) {
        var url_RunDetails = '/api/v1/runs/' + id + '/rundetails/';
        $.get(url_RunDetails, function (data) {
            //console.log(data);
            self.livedata.minIONname = data[0].minION_name;
            self.livedata.asicid = data[0].minKNOW_asic_id;
            self.livedata.scriptid = data[0].minKNOW_current_script;
            self.livedata.colours_string = data[0].minKNOW_colours_string;
            //self.livedata.computer=data[0].minKNOW_computer;
            //console.log(data[0]);
            document.getElementById('ComputerName').innerHTML = data[0].minKNOW_computer;
            document.getElementById('ScriptID').innerHTML = data[0].minKNOW_current_script;

        })
    }

    this.requestMessages = function () {
        var url_sincemessages = self.rundata.minION + 'messagessince/' + self.rundata.start_time + '/' + self.lasttime.toISOString() + "/";
        //console.log(url_sincemessages);
        $.get(url_sincemessages, function (data) {
            //console.log(data);
            stringtowrite = '<table class="table table-condensed"><tr><th>Message</th><th>Time</th></tr>';
            for (var i = 0; i < data.length; i++) {
                //stringtowrite=stringtowrite+'<div class="alert alert-info" role="alert">'+data[i].minKNOW_message + ' <p>(<i>' + new Date(data[i].minKNOW_message_timestamp) + '</i>) '+'</div>'
                stringtowrite = stringtowrite + '<tr><td>' + data[i].minKNOW_message + ' </td><td><i>' + new Date(data[i].minKNOW_message_timestamp) + '</i></td> ' + '</tr>'
            }
            stringtowrite = stringtowrite + '</table>';
            document.getElementById('Messages').innerHTML = stringtowrite;
        })
    }

    this.requestLiveRunStats = function (id) {
        //console.log('lastread ' + this.lastread);
        var url_livestats = '/api/v1/runs/' + id + '/runstats/' + this.lastread;
        $.get(url_livestats, function (data) {
            //console.log(data);
            if (data.length > 0) {
                self.needtoupdatecharts = true;
                self.lastread = data[data.length - 1].id;
                self.lasttime = new Date(data[data.length - 1].sample_time)
                //console.log(self.rundata.start_time);
                //console.log(self.lasttime.toISOString());
                self.requestMessages();
                for (var i = 0; i < data.length; i++) {
                    //console.log(data[i]);
                    timestamp = new Date(data[i].sample_time).getTime();
                    self.livedata.live_read_count = data[i].minKNOW_read_count;
                    self.livedata.voltage.push([timestamp, data[i].voltage_value]);
                    self.livedata.asictemp.push([timestamp, data[i].asic_temp]);
                    self.livedata.heatsinktemp.push([timestamp, data[i].heat_sink_temp]);
                    self.livedata.strand.push([timestamp, data[i].strand]);
                    self.livedata.good_single.push([timestamp, data[i].good_single]);
                    self.livedata.currpercentage = data[i].occupancy;
                    self.livedata.currstrand = data[i].strand;
                    self.livedata.percentage.push([timestamp, data[i].occupancy]);
                    self.livedata.yield_history.push([timestamp, data[i].event_yield]);
                    self.livedata.meanratio_history.push([timestamp, data[i].mean_ratio]);
                    self.livedata.instrand_history.push([timestamp, data[i].in_strand]);
                    self.livedata.openpore_history.push([timestamp, parseInt(data[i].open_pore)]);
                    var myStringArray = ["above", "adapter", "below", "good_single", "strand", "inrange", "multiple", "pending_mux_change", "saturated", "unavailable", "unblocking", "unclassified", "unknown"];
                    var arrayLength = myStringArray.length;
                    //console.log(parseInt(data[i][myStringArray[4]]));
                    for (var j = 0; j < arrayLength; j++) {
                        if (isNaN(data[i][myStringArray[j]])) {
                            self.livedata.pore_history[myStringArray[j]].push([timestamp, 0]);
                            //console.log("found a NAN");
                            //console.log(data[i][myStringArray[i]]);
                        } else {
                            self.livedata.pore_history[myStringArray[j]].push([timestamp, parseInt(data[i][myStringArray[j]])]);
                        }

                    }
                }
                if (self.needtoupdatecharts == true) {
                    self.updateLiveHistogram(data);
                    self.updateLiveYieldProjection();
                    self.updateLiveCumuYield();
                    self.updatePoreStats();
                }

            }
            //console.log(self.livedata);
        })
    };

    this.requestData = function () {

        var url_run = '/api/v1/runs/' + self.id;

        $.get(url_run, function (data) {
            //console.log(self.rundata);

            var barcodes = [];

            for (var i = 0; i < data.barcodes.length; i++) {
                barcodes.push(data.barcodes[i].name)
            }

            self.barcodes = barcodes.sort();
            self.rundata = data;
            self.updatetext(self.rundata);
            //console.log(self.rundata);
            self.updateBarcodeNavTab();

            //self.getlivedata(self.rundata);
        });

        self.requestRunDetails(self.id);
        self.requestLiveRunStats(self.id);
        self.requestSummaryByMinuteData(self.id);
        self.requestSummaryData(self.id);
        self.requestHistogramData(self.id);

        //console.log(self);

    };

}
