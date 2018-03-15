function updateReadsColumnBasedChart (chart, field, summaries, selectedBarcode) {

    var series = [];

    data = [];

    for (var readtype of Object.keys(summaries["All reads"])) {
        data.push(summaries["All reads"][readtype][field]['all']["data"][0]);
    }

    serie = {
        "name": "All reads",
        "data": data
    };

    series.push(serie);

    data = [];

    for (var readtype of Object.keys(summaries["All reads"])) {
        data.push(summaries["All reads"][readtype][field]['pass']["data"][0]);
    }

    serie = {
        "name": "All pass reads",
        "data": data
    };

    series.push(serie);

    if (field != "max_length") {

        data = [];

        for (var readtype of Object.keys(summaries["All reads"])) {
            data.push(summaries["All reads"][readtype][field]['fail']["data"][0]);
        }

        serie = {
            "name": "All fail reads",
            "data": data
        };

        series.push(serie);
    }

    // Include specific barcode if selected
    if (selectedBarcode !== "All reads") {

        data = [];
        //console.log(summaries[self.selectedBarcode]);
        for (var readtype of Object.keys(summaries[selectedBarcode])) {
            //console.log(readtype);
            //console.log("are we here?");
            data.push(summaries[selectedBarcode][readtype][field]['all']["data"][0]);
        }

        serie = {
            "name": selectedBarcode,
            "data": data
        };

        //console.log(serie);

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
}

function updateSummaryBasedCharts (summaries, selectedBarcode, charts) {

    for (var prop in charts) {

        updateReadsColumnBasedChart(charts[prop], prop, summaries, selectedBarcode);

    }

}

//function FlowcellPageApp() {
var FlowcellPageApp = {

    constructor () {

    },

    init: function () {

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
        this.makeHeatmapChart = makeHeatmapChart;
        this.makeStepLineChart = makeStepLineChart;

        this.write_run_data = write_run_data;

        this.lastread = 0;
        this.needtoupdatecharts = false;

        this.updatePoreChart = updatePoreChart;
        this.updateStepLineChart = updateStepLineChart;

        this.drawtaskbutton = drawtaskbutton;

        this.requestChannelSummaryData = requestChannelSummaryData;

        this.updateBarcodeNavTab = updateBarcodeNavTab;

        this.requestHistogramData = requestHistogramData;

        this.requestQualTime = requestQualTime;

        this.requestLengthTime = requestLengthTime;

        this.requestMaxLengthTime = requestMaxLengthTime;

        this.requestSpeed = requestSpeed;

        this.updateAssemblyCharts = updateAssemblyCharts;

        this.updateAssemblyBoxplot = updateAssemblyBoxplot;

        this.createAssemblyTable = createAssemblyTable;

        this.requestGfaData = requestGfaData;

        this.updatetext = updatetext;

        this.updateCumuBaseChart = updateCumuBaseChart;

        this.requestCumuBases = requestCumuBases;

        this.requestPafData = requestPafData;

        this.requestKraken = requestKraken;

        this.liveUpdateTasks = liveUpdateTasks;

        this.requestRunDetails = requestRunDetails;

        this.requestLiveRunStats = requestLiveRunStats;

        this.requestSummaryData = requestSummaryData;

        this.requestData = requestData;

        this.startTask = function (description, reference) {
            //console.log(description + " " + reference);
            var e = document.getElementById(description);
            if (e != null) {
                var strUser = e.options[e.selectedIndex].value;
            } else {
                var strUser = "null";
            }
            //console.log(strUser);
            $.ajaxSetup({
                beforeSend: function (xhr, settings) {
                    if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
                        xhr.setRequestHeader("X-CSRFToken", getCookie('csrftoken'));
                    }
                }
            });
            //var url_mask = "{% url 'set-task-detail-all' pk=12345 %}".replace(/12345/, reference.toString());
            var url_mask = "/api/v1/flowcells/12345/settask/".replace(/12345/, reference.toString());
            ;
            $.ajax({
                "type": "POST",
                "dataType": "json",
                "url": url_mask,
                "data": {
                    "job": description,
                    "reference": strUser,
                },
                "beforeSend": function (xhr, settings) {
                    //console.log("before send");
                    $.ajaxSettings.beforeSend(xhr, settings);
                },
                "success": function (result) {
                    //console.log(result);
                    $(".modal.in").modal("hide");
                    this.requestTasks(reference);

                },
                error: function (XMLHttpRequest, textStatus, errorThrown) {
                    //console.log("Status: " + textStatus);
                    //alert("Error: " + errorThrown);
                }

            })
        };

        this.chart_per_chrom_cov = this.makeChart(
            "per-chrom-cov",
            "Chromosome Coverage".toUpperCase(),
            "Chromosome Coverage".toUpperCase()
        );

        this.chart_per_chrom_avg = this.makeChart(
            "per-chrom-avg",
            "Read Length By Chromosome".toUpperCase(),
            "Read Length By Chromosome".toUpperCase()
        );

        this.chartChromosomeCoverage = this.makeStepLineChart(
            "chromosome-coverage",
            "Chromosome Coverage",
            "Coverage"
        );

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

        $('#chromosome-id-select').on('change', function () {
            this.updateStepLineChart(this.chartChromosomeCoverage, 0, 0);
        });

        $('#chromosome-coverage-left').on('click', function () {
            var min = this.chartChromosomeCoverage.xAxis[0].min;
            var max = this.chartChromosomeCoverage.xAxis[0].max;
            var delta = (max - min) / 2;

            this.updateStepLineChart(this.chartChromosomeCoverage, Math.round(min - delta, 0), Math.round(max - delta));
        });

        $('#chromosome-coverage-right').on('click', function () {
            var min = this.chartChromosomeCoverage.xAxis[0].min;
            var max = this.chartChromosomeCoverage.xAxis[0].max;
            var delta = (max - min) / 2;

            this.updateStepLineChart(this.chartChromosomeCoverage, Math.round(min + delta), Math.round(max + delta));
        });
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

        this.PoreShizzle = this.makeAreaPlot(
            'poreshizzle',
            'Pore States'.toUpperCase(),
            'Pore States'.toUpperCase()
        );

        console.log(">>> initializing monitorapp, requesting data");

        var inputFlowcellId = document.querySelector("#flowcell-id");

        this.flowcellId = inputFlowcellId.value;

        this.selectedBarcode = "All reads";

        this.requestData();

        //this.requestReference(this.id);

        //this.requestTasks(this.id);

        //setInterval(function () {
        //    this.requestData(this.id);
        //}.bind(this), 30000);

    }, // end of init

    //self: this;

    getSelectedBarcode: function () {

        if (!this.selectedBarcode) {
            return "All reads";
        } else {
            return this.selectedBarcode;
        }

    },


    /*
     * Each click on a barcode tab fires this function
     * that calls requestData and update all charts
     */
    updateChartsBasedOnBarcode: function (event) {
        self.selectedBarcode = event.target.innerText;
        self.requestData(self.id);
    },

    requestSummaryByMinuteData: function (id) {
        /*
         * Request summary by minute data
         */
        var url = "/api/v1/flowcells/" + id + "/summarybarcodebyminute";

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
    },

    updateSummaryByMinuteBasedCharts: function () {
        //self.updateAverageReadLengthOverTimeChart();
        //self.updateAverageQualityOverTimeChart();
        //self.updateCumulativeNumberOfReadsOverTimeChart();
        //self.updateCumulativeYieldOverTimeChart();
        //self.updateSequencingRateChart();
        //self.updateSequencingSpeedChart();
    },

    updateSequencingSpeedChart: function () {
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

        for (var i in self.rundata) {
            var starttime = new Date(Date.parse(self.rundata[i]['start_time']));
            var endtime = new Date(Date.parse(self.rundata[i]['last_read']));
            var name = self.rundata[i]['id']
            //chart.xAxis[0].addPlotBand({
            //    from: starttime,
            //    to: endtime,
            //    color: '#FCFFC5',
            //    id: name
            //});
            chart.xAxis[0].addPlotLine({
                value: starttime,
                color: 'black',
                dashStyle: 'dot',
                width: 2,
                //label: {
                //    text: name
                //}
            })
        }

    },

    updateSequencingRateChart: function () {
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

        for (var i in self.rundata) {
            var starttime = new Date(Date.parse(self.rundata[i]['start_time']));
            var endtime = new Date(Date.parse(self.rundata[i]['last_read']));
            var name = self.rundata[i]['id']
            //chart.xAxis[0].addPlotBand({
            //    from: starttime,
            //    to: endtime,
            //    color: '#FCFFC5',
            //    id: name
            //});
            chart.xAxis[0].addPlotLine({
                value: starttime,
                color: 'black',
                dashStyle: 'dot',
                width: 2,
                //label: {
                //    text: name
                //}
            })
        }

    },

    updateCumulativeYieldOverTimeChart: function () {
        var chart = self.chart_cumulative_yield_overtime;
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
                    summaries["All reads"][self.summaryByMinute[i].typename] = {};
                    summaries["All reads"][self.summaryByMinute[i].typename]['all'] = {
                        "lastCumulativeReadCount": 0,
                        "data": []
                    };
                    summaries["All reads"][self.summaryByMinute[i].typename]['pass'] = {
                        "lastCumulativeReadCount": 0,
                        "data": []
                    };
                    summaries["All reads"][self.summaryByMinute[i].typename]['fail'] = {
                        "lastCumulativeReadCount": 0,
                        "data": []
                    };
                }

                var cumulativeReadCount = self.summaryByMinute[i].total_length + summaries["All reads"][self.summaryByMinute[i].typename]['all'].lastCumulativeReadCount;
                var passcumulativeReadCount = self.summaryByMinute[i].pass_length + summaries["All reads"][self.summaryByMinute[i].typename]['pass'].lastCumulativeReadCount;
                var failcumulativeReadCount = self.summaryByMinute[i].total_length - self.summaryByMinute[i].pass_length + summaries["All reads"][self.summaryByMinute[i].typename]['fail'].lastCumulativeReadCount;
                var sample_time = new Date(self.summaryByMinute[i].sample_time);// - self.flowcellstart;

                var point = {
                    x: sample_time,
                    y: cumulativeReadCount
                }

                var passpoint = {
                    x: sample_time,
                    y: passcumulativeReadCount
                }

                var failpoint = {
                    x: sample_time,
                    y: failcumulativeReadCount
                }
                summaries["All reads"][self.summaryByMinute[i].typename]['all'].lastCumulativeReadCount = cumulativeReadCount;
                summaries["All reads"][self.summaryByMinute[i].typename]['all'].data.push(point);
                summaries["All reads"][self.summaryByMinute[i].typename]['pass'].lastCumulativeReadCount = passcumulativeReadCount;
                summaries["All reads"][self.summaryByMinute[i].typename]['pass'].data.push(passpoint);
                summaries["All reads"][self.summaryByMinute[i].typename]['fail'].lastCumulativeReadCount = failcumulativeReadCount;
                summaries["All reads"][self.summaryByMinute[i].typename]['fail'].data.push(failpoint);

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
                for (var qual in summaries[barcode][readtype]) {
                    chart.addSeries({
                        name: barcode + " - " + readtype + " - " + qual,
                        data: summaries[barcode][readtype][qual]["data"]
                    });
                }
            }
        }

        for (var i in self.rundata) {
            var starttime = new Date(Date.parse(self.rundata[i]['start_time']));// - self.flowcellstart;
            var endtime = new Date(Date.parse(self.rundata[i]['last_read']));//- self.flowcellstart;
            var name = self.rundata[i]['id']
            //chart.xAxis[0].addPlotBand({
            //    from: starttime,
            //    to: endtime,
            //    color: '#FCFFC5',
            //    id: name
            //});
            chart.xAxis[0].addPlotLine({
                value: starttime,
                color: 'black',
                dashStyle: 'dot',
                width: 2,
                //label: {
                //    text: name
                //}
            })
        }

    },

    updateCumulativeNumberOfReadsOverTimeChart: function () {
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
                    summaries["All reads"][self.summaryByMinute[i].typename] = {};
                    summaries["All reads"][self.summaryByMinute[i].typename]['all'] = {
                        "lastCumulativeReadCount": 0,
                        "data": []
                    };
                    summaries["All reads"][self.summaryByMinute[i].typename]['pass'] = {
                        "lastCumulativeReadCount": 0,
                        "data": []
                    };
                    summaries["All reads"][self.summaryByMinute[i].typename]['fail'] = {
                        "lastCumulativeReadCount": 0,
                        "data": []
                    };
                }

                var cumulativeReadCount = self.summaryByMinute[i].read_count + summaries["All reads"][self.summaryByMinute[i].typename]['all'].lastCumulativeReadCount;
                var passcumulativeReadCount = self.summaryByMinute[i].pass_count + summaries["All reads"][self.summaryByMinute[i].typename]['pass'].lastCumulativeReadCount;
                var failcumulativeReadCount = self.summaryByMinute[i].read_count - self.summaryByMinute[i].pass_count + summaries["All reads"][self.summaryByMinute[i].typename]['fail'].lastCumulativeReadCount;
                var sample_time = new Date(self.summaryByMinute[i].sample_time);// - self.flowcellstart;

                var point = {
                    x: sample_time,
                    y: cumulativeReadCount
                }

                var passpoint = {
                    x: sample_time,
                    y: passcumulativeReadCount
                }

                var failpoint = {
                    x: sample_time,
                    y: failcumulativeReadCount
                }
                summaries["All reads"][self.summaryByMinute[i].typename]['all'].lastCumulativeReadCount = cumulativeReadCount;
                summaries["All reads"][self.summaryByMinute[i].typename]['all'].data.push(point);
                summaries["All reads"][self.summaryByMinute[i].typename]['pass'].lastCumulativeReadCount = passcumulativeReadCount;
                summaries["All reads"][self.summaryByMinute[i].typename]['pass'].data.push(passpoint);
                summaries["All reads"][self.summaryByMinute[i].typename]['fail'].lastCumulativeReadCount = failcumulativeReadCount;
                summaries["All reads"][self.summaryByMinute[i].typename]['fail'].data.push(failpoint);

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
                for (var qual in summaries[barcode][readtype]) {
                    chart.addSeries({
                        name: barcode + " - " + readtype + " - " + qual,
                        data: summaries[barcode][readtype][qual]["data"]
                    });
                }
            }
        }

        for (var i in self.rundata) {
            var starttime = new Date(Date.parse(self.rundata[i]['start_time']));// - self.flowcellstart;
            var endtime = new Date(Date.parse(self.rundata[i]['last_read']));//- self.flowcellstart;
            var name = self.rundata[i]['id']
            //chart.xAxis[0].addPlotBand({
            //    from: starttime,
            //    to: endtime,
            //    color: '#FCFFC5',
            //    id: name
            //});
            chart.xAxis[0].addPlotLine({
                value: starttime,
                color: 'black',
                dashStyle: 'dot',
                width: 2,
                //label: {
                //    text: name
                //}
            })
        }

    },

    updateAverageQualityOverTimeChart: function () {

        var chart = self.average_quality_overtime;
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

            var average_quality = self.summaryByMinute[i].quality_sum / self.summaryByMinute[i].read_count;
            var pass_average_quality = self.summaryByMinute[i].pass_quality_sum / self.summaryByMinute[i].pass_count;
            var fail_average_quality = (self.summaryByMinute[i].quality_sum - self.summaryByMinute[i].pass_quality_sum) / (self.summaryByMinute[i].read_count - self.summaryByMinute[i].pass_count);

            var sample_time = new Date(self.summaryByMinute[i].sample_time);

            var point = {
                x: sample_time,
                y: average_quality
            };

            var pass_point = {
                x: sample_time,
                y: pass_average_quality
            };

            var fail_point = {
                x: sample_time,
                y: fail_average_quality
            };

            if (self.summaryByMinute[i].barcodename === "All reads") {
                if (summaries["All reads"][self.summaryByMinute[i].typename] === undefined) {
                    summaries["All reads"][self.summaryByMinute[i].typename] = [];
                    summaries["All reads"][self.summaryByMinute[i].typename]['all'] = [];
                    summaries["All reads"][self.summaryByMinute[i].typename]['pass'] = [];
                    summaries["All reads"][self.summaryByMinute[i].typename]['fail'] = [];
                }

                summaries["All reads"][self.summaryByMinute[i].typename]['all'].push(point);
                summaries["All reads"][self.summaryByMinute[i].typename]['pass'].push(pass_point);
                summaries["All reads"][self.summaryByMinute[i].typename]['fail'].push(fail_point);
            }

            if (self.summaryByMinute[i].barcodename === selectedBarcode && selectedBarcode !== "All reads") {
                if (summaries[selectedBarcode][self.summaryByMinute[i].typename] === undefined) {
                    summaries[selectedBarcode][self.summaryByMinute[i].typename] = [];
                    summaries[selectedBarcode][self.summaryByMinute[i].typename]['all'] = [];
                    summaries[selectedBarcode][self.summaryByMinute[i].typename]['pass'] = [];
                    summaries[selectedBarcode][self.summaryByMinute[i].typename]['fail'] = [];


                }

                summaries[selectedBarcode][self.summaryByMinute[i].typename]['all'].push(point);
                summaries[selectedBarcode][self.summaryByMinute[i].typename]['pass'].push(pass_point);
                summaries[selectedBarcode][self.summaryByMinute[i].typename]['fail'].push(fail_point);
            }

        }

        for (var barcode in summaries) {
            for (var readtype in summaries[barcode]) {
                for (var qual in summaries[barcode][readtype]) {
                    chart.addSeries({
                        name: barcode + " - " + readtype + " - " + qual,
                        data: summaries[barcode][readtype][qual]
                    });
                }
            }
        }

        for (var i in self.rundata) {
            var starttime = new Date(Date.parse(self.rundata[i]['start_time']));
            var endtime = new Date(Date.parse(self.rundata[i]['last_read']));
            var name = self.rundata[i]['id']
            //chart.xAxis[0].addPlotBand({
            //    from: starttime,
            //    to: endtime,
            //    color: '#FCFFC5',
            //    id: name
            //});
            chart.xAxis[0].addPlotLine({
                value: starttime,
                color: 'black',
                dashStyle: 'dot',
                width: 2,
                //label: {
                //    text: name
                //}
            })
        }
    },

    updateAverageReadLengthOverTimeChart: function () {

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

        for (var i in self.rundata) {
            var starttime = new Date(Date.parse(self.rundata[i]['start_time']));
            var endtime = new Date(Date.parse(self.rundata[i]['last_read']));
            var name = self.rundata[i]['id']
            //chart.xAxis[0].addPlotBand({
            //    from: starttime,
            //    to: endtime,
            //    color: '#FCFFC5',
            //    id: name
            //});
            chart.xAxis[0].addPlotLine({
                value: starttime,
                color: 'black',
                dashStyle: 'dot',
                width: 2,
                //label: {
                //    text: name
                //}
            })
        }
    },


    updatePoreStats: function () {
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
    },

    parseporehist: function (descriptions, counts) {
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
    },

    updateLiveCumuYield: function () {
        //console.log(self.LiveCumuYield);
        if (self.LiveCumuYield.series.length < 1) {
            self.LiveCumuYield.addSeries({
                data: self.livedata.yield_history
            });
            self.LiveCumuYield.series[0].update({name: "Events"}, false);
        } else {
            self.LiveCumuYield.series[0].setData(self.livedata.yield_history);
            self.LiveCumuYield.series[0].update({name: "Events"}, false);
        }
        self.LiveCumuYield.redraw();
        self.LiveCumuYield.reflow();
        if (self.LiveOccupancy.series.length < 1) {
            self.LiveOccupancy.addSeries({
                data: self.livedata.percentage
            });
            self.LiveOccupancy.series[0].update({name: "% Occupancy"}, false);
        } else {
            self.LiveOccupancy.series[0].setData(self.livedata.percentage);
            self.LiveOccupancy.series[0].update({name: "% Occupancy"}, false);
        }
        self.LiveOccupancy.redraw();
        self.LiveOccupancy.reflow();
        if (self.LiveInStrand.series.length < 2) {
            self.LiveInStrand.addSeries({data: self.livedata.strand});
            self.LiveInStrand.series[0].update({name: "In Strand"}, false);
            self.LiveInStrand.addSeries({data: self.livedata.good_single});
            self.LiveInStrand.series[1].update({name: "Single Pore"}, false);

        } else {
            self.LiveInStrand.series[0].setData(self.livedata.strand);
            self.LiveInStrand.series[0].update({name: "In Strand"}, false);
            self.LiveInStrand.series[1].setData(self.livedata.good_single);
            self.LiveInStrand.series[1].update({name: "Single Pore"}, false);
        }
        self.LiveInStrand.redraw();
        self.LiveInStrand.reflow();
        if (self.LivePoreState.series.length < 1) {
            self.LivePoreState.addSeries({data: self.livedata.instrand_history});
            self.LivePoreState.addSeries({data: self.livedata.openpore_history});
        } else {
            self.LivePoreState.series[0].setData(self.livedata.instrand_history);
            self.LivePoreState.series[1].setData(self.livedata.openpore_history);
        }
        self.LivePoreState.series[0].update({name: "In Strand"}, false);
        self.LivePoreState.series[1].update({name: "Open Pore"}, false);
        self.LivePoreState.redraw();
        self.LivePoreState.reflow();
        if (self.LiveCurrentRatio.series.length < 1) {
            self.LiveCurrentRatio.addSeries({data: self.livedata.meanratio_history});
        } else {
            self.LiveCurrentRatio.series[0].setData(self.livedata.meanratio_history);
        }
        self.LiveCurrentRatio.series[0].update({name: "Current Ratio"}, false);
        self.LiveCurrentRatio.redraw();
        self.LiveCurrentRatio.reflow();
        if (self.LiveTemperature.series.length < 1) {
            self.LiveTemperature.addSeries({data: self.livedata.asictemp});
            self.LiveTemperature.addSeries({data: self.livedata.heatsinktemp});
        } else {
            self.LiveTemperature.series[0].setData(self.livedata.asictemp);
            self.LiveTemperature.series[1].setData(self.livedata.heatsinktemp);
        }
        self.LiveTemperature.series[0].update({name: "Asic Temp"}, false);
        self.LiveTemperature.series[1].update({name: "HeatSink Temp"}, false);
        self.LiveTemperature.redraw();
        self.LiveTemperature.reflow();
        if (self.LiveVoltage.series.length < 1) {
            self.LiveVoltage.addSeries({data: self.livedata.voltage});
        } else {
            self.LiveVoltage.series[0].setData(self.livedata.voltage);
        }
        self.LiveVoltage.series[0].update({name: "Voltage"}, false);
        self.LiveVoltage.redraw();
        self.LiveVoltage.reflow();

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
    },

    converttobases: function (data, seqspeed) {
        if (Number(self.livedata.scalingfactor) > Number(0)) {
            //console.log("returning scaling factor" + self.livedata.scalingfactor);
            scaling = Number(self.livedata.scalingfactor);
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

    requestMessages: function () {

        console.log(">>> antes do rundetails");
        if (!self.rundetails) {
            return;
        }

        var url_sincemessages = self.rundetails[0]["minION"] + 'messagessince/' + self.rundetails[0]["minKNOW_start_time"] + '/' + self.lasttime.toISOString() + "/";

        console.log(url_sincemessages);

        $.get(url_sincemessages, function (data) {
            console.log(data);
            stringtowrite = '<table class="table table-condensed"><tr><th>Message</th><th>Time</th></tr>';
            for (var i = 0; i < data.length; i++) {
                //stringtowrite=stringtowrite+'<div class="alert alert-info" role="alert">'+data[i].minKNOW_message + ' <p>(<i>' + new Date(data[i].minKNOW_message_timestamp) + '</i>) '+'</div>'
                stringtowrite = stringtowrite + '<tr><td>' + data[i].minKNOW_message + ' </td><td><i>' + new Date(data[i].minKNOW_message_timestamp) + '</i></td> ' + '</tr>'
            }
            stringtowrite = stringtowrite + '</table>';
            document.getElementById('Messages').innerHTML = stringtowrite;
        })
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

    },

    updateLiveHistogram: function (data) {
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
    },


    calculatereadtoeventscaling: function () {
        var totalyield = 0;
        var readcount = 0;
        if (self.summary !== null) {
            //if ('All reads' in self.summary) {
            for (var readtype in self.summary["All reads"]) {
                //console.log(self.summary['All reads'][readtype]);
                totalyield = totalyield + parseInt(self.summary['All reads'][readtype]['yield']['data']);
                readcount = readcount + parseInt(self.summary['All reads'][readtype]['read_count']['data']);
            }
            //console.log("yieldhistory length");
            //console.log("test" + self.livedata.live_read_count);
            if (self.livedata.yield_history.length > 1) {
                self.livedata.scalingfactor = (totalyield / readcount) / (self.livedata.yield_history[self.livedata.yield_history.length - 1][1] / self.livedata.live_read_count);
            }
        }
    },

    requestReference: function (id) {
        var url = "/api/v1/reference/";
        $.get(url, function (data) {
            var references = [];
            for (var i = 0; i < data.length; i++) {
                //console.log(data[i]);
                references.push(data[i]);
            }
            //console.log(references);
            self.references = references;
        });
    },

    requestTasks: function (id) {
        var url = "/api/v1/flowcells/" + id + "/tasks/";
        $.get(url, function (data) {
            //console.log(data);
            var tasks = [];
            for (var i = 0; i < data.length; i++) {
                tasks.push(data[i]);
            }
            self.tasks = tasks;
            self.updateTasks(id); //really this only needs to run once!
        });
    },


    updateTasks: function (id) {
        //console.log("####id is" + id);
        var taskstring = "";
        for (var i = 0; i < self.tasks.length; i++) {
            //console.log(self.tasks[i]);
            if (self.tasks[i].hasOwnProperty("job_details")) {
                colour = 'bg-green';
                message = 'Reads processed:' + self.tasks[i]["job_details"]["read_count"];
                percentage = 50;
                message2 = 'X% of uploaded reads are processed';
                icon = 'fa fa-refresh fa-spin fa-fw';
                running = true;

            } else {
                colour = 'bg-light-blue';
                message = 'Task Not Running.';
                percentage = 0;
                message2 = "Click to start a " + self.tasks[i]["description"] + ' task.';
                icon = 'fa fa-refresh fa-fw';
                running = false;
            }
            taskstring = this.drawtaskbutton(taskstring, colour, icon, self.tasks[i]["description"], message, percentage, message2, i, self.tasks[i]["long_description"], self.tasks[i]["reference"], self.tasks[i]["name"], self.tasks[i]["transcriptome"]);
        }
        ;
        document.getElementById('tasks').innerHTML = taskstring;
        for (var i = 0; i < self.tasks.length; i++) {
            var buttonname = "#button" + self.tasks[i]["name"];
            //console.log("Button name to look for:" + buttonname);
            $(buttonname).click(function (e) {
                var idClicked = e.target.id;
                //console.log(idClicked);
                self.startTask(idClicked.substring(6), id);
            });
        }
        ;
        self.liveUpdateTasks(id);

    },


    //this.requestMappedChromosomes = requestMappedChromosomes;

    updateCoverageBasedCharts: function (chart, field) {
        var summarycoverage = this.summarycoverage;
        var series = [];
        var categories = [];
        for (var barcode of Object.keys(summarycoverage)) {
    //            console.log(summarycoverage[barcode]);
            data = [];
            for (var readtype of Object.keys(summarycoverage[barcode])) {

                for (var chromosome of Object.keys(summarycoverage[barcode][readtype])) {
                    categories.push(chromosome);
    //                    console.log(summarycoverage[barcode][readtype][chromosome]['coverage']["data"]);
                    data.push(summarycoverage[barcode][readtype][chromosome][field]["data"]);
                }
            }
            serie = {
                "name": barcode + ' ' + readtype,
                "data": data
            };
            series.push(serie);
        }
        ;
        chart.xAxis[0].setCategories(categories);
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
    }

}
