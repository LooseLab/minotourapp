function MonitorAPP() {
    var livedata = new Array();


    this.init = function () {
        this.requestData();

        setInterval(function () {
            this.requestData();
        }.bind(this), 30000);
    };

    this.updatecounters = function (live) {
        //var thing = document.getElementById("livenum");
        //thing.innerHTML = live;
        //var thing2 = document.getElementById("livenumruns");
        //thing2.innerHTML = "You have " + live + " live runs.";
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

function MinotourFlowCellApp() {

    this.chart_reads_called = null;
    this.chart_yield = null;
    this.chart_average_read_length = null;
    this.chart_maximum_read_length = null;
    //this.average_read_lengths_overtime = null;
    this.average_read_lengths_overtime_new = null;
    this.xy_scat_length = null;
    this.trans_top100 = null;
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

    this.chart_per_chrom_cov = null;

    this.coveragedata = new Array();
    this.coveragedata.read_type = new Array();

    this.makeChart = makeChart;
    this.makeChart2 = makeChart2;
    this.makeChart3 = makeChart3;
    this.makeChart4 = makeChart4;
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
                self.requestTasks(reference);

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

    /*this.average_read_lengths_overtime = this.makeChart2(
        "average-read-lengths-overtime",
        "average read length over time".toUpperCase(),
        "average read length".toUpperCase()
    );*/

    this.average_read_lengths_overtime_new = this.makeChart2(
        "average-read-lengths-overtime-new",
        "average read length over time".toUpperCase(),
        "average read length".toUpperCase()
    );

    /*this.average_quality_overtime = this.makeChart2(
        "average-quality-overtime",
        "average quality over time".toUpperCase(),
        "average read quality score".toUpperCase()
    );*/

    this.average_quality_overtime_new = this.makeChart2(
        "average-quality-overtime-new",
        "average quality over time".toUpperCase(),
        "average read quality score".toUpperCase()
    );

    /*this.chart_cumulative_number_reads_overtime = this.makeChart2(
        "cumulative-number-reads-overtime",
        "cumulative reads".toUpperCase(),
        "cumulative reads".toUpperCase()
    );*/

    this.chart_cumulative_number_reads_overtime_new = this.makeChart2(
        "cumulative-number-reads-overtime-new",
        "cumulative reads".toUpperCase(),
        "cumulative reads".toUpperCase()
    );

    /*
    this.chart_cumulative_yield_overtime = this.makeChart2(
        "cumulative-yield-overtime",
        "cumulative bases".toUpperCase(),
        "cumulative bases".toUpperCase()
    );
    */

    this.chart_cumulative_yield_overtime_new = this.makeChart2(
        "cumulative-yield-overtime-new",
        "cumulative bases".toUpperCase(),
        "cumulative bases".toUpperCase()
    );

    this.chartSequencingRate_new = this.makeChart2(
        "sequencing-rate-new",
        "sequencing rate".toUpperCase(),
        "bases/second".toUpperCase()
    );

    this.chartSequencingSpeed_new = this.makeChart2(
        "sequencing-speed-new",
        "sequencing speed".toUpperCase(),
        "bases/channel/second".toUpperCase()
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

    this.chartReadsPerPore = this.makeHeatmapChart(
        "reads-per-pore",
        "Reads per Channel".toUpperCase(),
        ""
    );

    this.chartBasesPerPore = this.makeHeatmapChart(
        "bases-per-pore",
        "bases (kb) per Channel".toUpperCase(),
        ""
    );

    this.chartChromosomeCoverage = this.makeStepLineChart(
        "chromosome-coverage",
        "Chromosome Coverage",
        "Coverage"
    );

    $('#chromosome-id-select').on('change', function () {
        self.updateStepLineChart(self.chartChromosomeCoverage, 0, 0);
    });

    $('#chromosome-coverage-left').on('click', function () {
        var min = self.chartChromosomeCoverage.xAxis[0].min;
        var max = self.chartChromosomeCoverage.xAxis[0].max;
        var delta = (max - min) / 2;

        self.updateStepLineChart(self.chartChromosomeCoverage, Math.round(min - delta, 0), Math.round(max - delta));
    });

    $('#chromosome-coverage-right').on('click', function () {
        var min = self.chartChromosomeCoverage.xAxis[0].min;
        var max = self.chartChromosomeCoverage.xAxis[0].max;
        var delta = (max - min) / 2;

        self.updateStepLineChart(self.chartChromosomeCoverage, Math.round(min + delta), Math.round(max + delta));
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


    var self = this;

    this.init = function () {
        self.id = document.getElementById("flowcell-id").innerText;
        self.selectedBarcode = "All reads";

        this.requestData();

        this.requestReference(self.id);
        this.requestTasks(self.id);

        setInterval(function () {
            this.requestData();
        }.bind(this), 30000);

    };

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

        var sortedBarcodes = self.barcodes;

        for (var i = 0; i < sortedBarcodes.length; i++) {
            var li = document.createElement("li");
            var a = document.createElement("a");
            a.onclick = self.updateChartsBasedOnBarcode;
            a.href = "#";
            a.text = sortedBarcodes[i];

            if (sortedBarcodes[i] === self.selectedBarcode) {
                li.classList.add("active");
            }

            li.appendChild(a);
            ul.appendChild(li);
        }
    };

    this.requestCumuBases = function (id){
        var url = "/api/v1/flowcells/" + id + "/summarybarcodebyminute_bases";
        $.get(url,function (data){
            self.cumulativereadsbybarcode=data;
            self.updateCumuBaseChart();
        })

    };

    this.requestSpeed = function (id){
        var url = "/api/v1/flowcells/" + id + "/summarybarcodebyminute_speed";
        $.get(url,function (data){
            self.speedreadsbybarcode=data;
            self.updateSpeedBaseChart();
        })

    };

    this.requestQualTime = function (id){
        var url = "/api/v1/flowcells/" + id + "/summarybarcodebyminute_quality";
        $.get(url,function (data){
            self.qualityovertime=data;
            self.updateQualChart();
        })
    };

    this.requestLengthTime = function (id){
        var url = "/api/v1/flowcells/" + id + "/summarybarcodebyminute_length";
        $.get(url,function (data){
            self.lengthovertime=data;
            self.updateLengthChart();
        })
    };

    this.updateLengthChart = function () {
        var chart = self.average_read_lengths_overtime_new;
        var selectedBarcode = self.selectedBarcode;
        // Remove previous series
        while (chart.series.length > 0) {
            chart.series[0].remove();
        }
        for (var typeName of Object.keys(self.lengthovertime[selectedBarcode])) {
            for (var status of Object.keys(self.lengthovertime[selectedBarcode][typeName])) {
                chart.addSeries({
                    name: selectedBarcode + " - " + typeName + " - " + status,
                    data: self.lengthovertime[selectedBarcode][typeName][status]
                });
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
    };

    this.updateQualChart = function () {
        var chart = self.average_quality_overtime_new;
        var selectedBarcode = self.selectedBarcode;
        // Remove previous series
        while (chart.series.length > 0) {
            chart.series[0].remove();
        }
        for (var typeName of Object.keys(self.qualityovertime[selectedBarcode])) {
            for (var status of Object.keys(self.qualityovertime[selectedBarcode][typeName])) {
                chart.addSeries({
                    name: selectedBarcode + " - " + typeName + " - " + status,
                    data: self.qualityovertime[selectedBarcode][typeName][status]
                });
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
    };

    this.updateSpeedBaseChart = function () {
        var chart = self.chartSequencingRate_new;
        var chart2 = self.chartSequencingSpeed_new;
        var selectedBarcode = self.selectedBarcode;


        // Remove previous series
        while (chart.series.length > 0) {
            chart.series[0].remove();
        }

        // Remove previous series
        while (chart2.series.length > 0) {
            chart2.series[0].remove();
        }

        for (var typeName of Object.keys(self.speedreadsbybarcode[selectedBarcode]["speed"])) {
            for (var status of Object.keys(self.speedreadsbybarcode[selectedBarcode]["speed"][typeName])) {
                chart2.addSeries({
                    name: selectedBarcode + " - " + typeName + " - " + status,
                    data: self.speedreadsbybarcode[selectedBarcode]["speed"][typeName][status]
                });
            }
        }

        for (var typeName of Object.keys(self.speedreadsbybarcode[selectedBarcode]["rate"])) {
            for (var status of Object.keys(self.speedreadsbybarcode[selectedBarcode]["rate"][typeName])) {
                chart.addSeries({
                    name: selectedBarcode + " - " + typeName + " - " + status,
                    data: self.speedreadsbybarcode[selectedBarcode]["rate"][typeName][status]
                });
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
            chart2.xAxis[0].addPlotLine({
                value: starttime,
                color: 'black',
                dashStyle: 'dot',
                width: 2,
                //label: {
                //    text: name
                //}
            })
        }
    };


    this.updateCumuBaseChart = function () {
        var chart = self.chart_cumulative_yield_overtime_new;
        var chart2 = self.chart_cumulative_number_reads_overtime_new;
        var selectedBarcode = self.selectedBarcode;


        // Remove previous series
        while (chart.series.length > 0) {
            chart.series[0].remove();
        }

        // Remove previous series
        while (chart2.series.length > 0) {
            chart2.series[0].remove();
        }

        for (var typeName of Object.keys(self.cumulativereadsbybarcode[selectedBarcode]["bases"])) {
            for (var status of Object.keys(self.cumulativereadsbybarcode[selectedBarcode]["bases"][typeName])) {
                chart.addSeries({
                    name: selectedBarcode + " - " + typeName + " - " + status,
                    data: self.cumulativereadsbybarcode[selectedBarcode]["bases"][typeName][status]
                });
            }
        }

        for (var typeName of Object.keys(self.cumulativereadsbybarcode[selectedBarcode]["reads"])) {
            for (var status of Object.keys(self.cumulativereadsbybarcode[selectedBarcode]["reads"][typeName])) {
                chart2.addSeries({
                    name: selectedBarcode + " - " + typeName + " - " + status,
                    data: self.cumulativereadsbybarcode[selectedBarcode]["reads"][typeName][status]
                });
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
            chart2.xAxis[0].addPlotLine({
                value: starttime,
                color: 'black',
                dashStyle: 'dot',
                width: 2,
                //label: {
                //    text: name
                //}
            })
        }
    };


    this.requestSummaryByMinuteData = function (id) {
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
    };

    this.updateSummaryByMinuteBasedCharts = function () {
        //self.updateAverageReadLengthOverTimeChart();
        //self.updateAverageQualityOverTimeChart();
        //self.updateCumulativeNumberOfReadsOverTimeChart();
        //self.updateCumulativeYieldOverTimeChart();
        //self.updateSequencingRateChart();
        //self.updateSequencingSpeedChart();
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

    };

    this.updateCumulativeYieldOverTimeChart = function () {
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

    };

    this.updateAverageQualityOverTimeChart = function () {

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
    };

    this.requestSummaryData = function (id) {
        /*
         * Request summary by barcode data
         */
        var url = "/api/v1/flowcells/" + id + "/summarybarcode";

        $.get(url, function (data) {
            if (data.length > 0) {
                var summary = {};
                for (var i = 0; i < data.length; i++) {
                    var item = data[i];

                    if (summary[item.barcodename] === undefined) {
                        summary[item.barcodename] = {};
                    }

                    if (summary[item.barcodename][item.typename] === undefined) {
                        summary[item.barcodename][item.typename] = {};
                        summary[item.barcodename][item.typename]["read_count"] = 0;
                        summary[item.barcodename][item.typename]["pass_count"] = 0;
                        summary[item.barcodename][item.typename]["yield"] = 0;
                        summary[item.barcodename][item.typename]["pass_yield"] = 0;
                        summary[item.barcodename][item.typename]["max_length"] = 0;
                        summary[item.barcodename][item.typename]["pass_max_length"] = 0;

                    }
                    summary[item.barcodename][item.typename]["read_count"] += item.read_count;
                    summary[item.barcodename][item.typename]["pass_count"] += item.pass_count;
                    summary[item.barcodename][item.typename]["yield"] += item.total_length;
                    summary[item.barcodename][item.typename]["pass_yield"] += item.pass_length;
                    if (item.max_length > summary[item.barcodename][item.typename]["max_length"]) {
                        summary[item.barcodename][item.typename]["max_length"] = item.max_length;
                    }
                    if (item.pass_max_length > summary[item.barcodename][item.typename]["pass_max_length"]) {
                        summary[item.barcodename][item.typename]["pass_max_length"] = item.pass_max_length;
                    }


                }
                var summaries = {};
                for (var barcodename in summary) {
                    if (summaries[barcodename] === undefined) {
                        summaries[barcodename] = {};
                    }
                    for (var typename in summary[barcodename]) {
                        if (summaries[barcodename][typename] === undefined) {
                            summaries[barcodename][typename] = {
                                "read_count": null,
                                "yield": null,
                                "average_read_length": null,
                                "max_length": null
                            };

                            summaries[barcodename][typename]["read_count"] = {};
                            summaries[barcodename][typename]["read_count"]['all'] = {
                                "name": typename,
                                "data": [summary[barcodename][typename]["read_count"]],//[summaries[item.barcodename][item.typename]["read_count"]["data"] + item.read_count],
                                "animation": false
                            };
                            summaries[barcodename][typename]["read_count"]['pass'] = {
                                "name": typename,
                                "data": [summary[barcodename][typename]["pass_count"]],//[summaries[item.barcodename][item.typename]["read_count"]["data"] + item.read_count],
                                "animation": false
                            };
                            summaries[barcodename][typename]["read_count"]['fail'] = {
                                "name": typename,
                                "data": [(summary[barcodename][typename]["read_count"] - summary[barcodename][typename]["pass_count"])],//[summaries[item.barcodename][item.typename]["read_count"]["data"] + item.read_count],
                                "animation": false
                            };

                            summaries[barcodename][typename]["yield"] = {};
                            summaries[barcodename][typename]["yield"]['all'] = {
                                "name": typename,
                                "data": [summary[barcodename][typename]["yield"]],
                                "animation": false
                            };
                            summaries[barcodename][typename]["yield"]['pass'] = {
                                "name": typename,
                                "data": [summary[barcodename][typename]["pass_yield"]],
                                "animation": false
                            };
                            summaries[barcodename][typename]["yield"]['fail'] = {
                                "name": typename,
                                "data": [(summary[barcodename][typename]["yield"] - summary[barcodename][typename]["pass_yield"])],
                                "animation": false
                            };

                            summaries[barcodename][typename]["average_read_length"] = {};
                            summaries[barcodename][typename]["average_read_length"]['all'] = {
                                "name": typename,
                                "data": [summary[barcodename][typename]["yield"] / summary[barcodename][typename]["read_count"]],
                                "animation": false
                            };
                            summaries[barcodename][typename]["average_read_length"]['pass'] = {
                                "name": typename,
                                "data": [summary[barcodename][typename]["pass_yield"] / summary[barcodename][typename]["pass_count"]],
                                "animation": false
                            };
                            summaries[barcodename][typename]["average_read_length"]['fail'] = {
                                "name": typename,
                                "data": [(summary[barcodename][typename]["yield"] - summary[barcodename][typename]["pass_yield"]) / (summary[barcodename][typename]["read_count"] - summary[barcodename][typename]["pass_count"])],
                                "animation": false
                            };

                            summaries[barcodename][typename]["max_length"] = {};
                            summaries[barcodename][typename]["max_length"]['all'] = {
                                "name": typename,
                                "data": [summary[barcodename][typename]["max_length"]],
                                "animation": false
                            };
                            summaries[barcodename][typename]["max_length"]['pass'] = {
                                "name": typename,
                                "data": [summary[barcodename][typename]["pass_max_length"]],
                                "animation": false
                            };
                            summaries[barcodename][typename]["max_length"]['fail'] = {
                                "name": typename,
                                "data": [summary[barcodename][typename]["max_length"]],
                                "animation": false
                            };

                        }
                    }

                }

                self.summary = summaries;

                self.updateSummaryBasedCharts();
                //console.log(self.livedata);

            }
        });

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

    this.updateReadsColumnBasedChart = function (chart, field) {
        var summaries = self.summary;

        var series = [];
        //console.log(summaries);
        // Always include all reads
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
        if (self.selectedBarcode !== "All reads") {

            data = [];
            //console.log(summaries[self.selectedBarcode]);
            for (var readtype of Object.keys(summaries[self.selectedBarcode])) {
                //console.log(readtype);
                //console.log("are we here?");
                data.push(summaries[self.selectedBarcode][readtype][field]['all']["data"][0]);
            }

            serie = {
                "name": self.selectedBarcode,
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
    };

    this.requestHistogramData = function (id) {
        /*
         * Request histogram data
         */

        var url = "/api/v1/flowcells/" + id + "/histogramsummary";

        $.get(url, function (data) {
            //console.log(data);
            if (data.length > 0) {

                var ordered_data = data.sort(function (a, b) {
                    return a.bin_width - b.bin_width;
                });

                var summary = {};
                //console.log(ordered_data);
                for (var i = 0; i < ordered_data.length; i++) {

                    var item = ordered_data[i];
                    //console.log(item);
                    if (summary[item.barcode_name] === undefined) {
                        summary[item.barcode_name] = {};
                    }

                    if (summary[item.barcode_name][item.read_type_name] === undefined) {
                        summary[item.barcode_name][item.read_type_name] = {};
                    }

                    if (summary[item.barcode_name][item.read_type_name][parseInt(item.bin_width)] === undefined) {
                        summary[item.barcode_name][item.read_type_name][parseInt(item.bin_width)] = {};
                        summary[item.barcode_name][item.read_type_name][parseInt(item.bin_width)]['read_count'] = 0;
                        summary[item.barcode_name][item.read_type_name][parseInt(item.bin_width)]['read_length'] = 0;
                    }
                    summary[item.barcode_name][item.read_type_name][parseInt(item.bin_width)]['read_count'] += item.read_count;
                    summary[item.barcode_name][item.read_type_name][parseInt(item.bin_width)]['read_length'] += item.read_length;
                }
                ;
                //console.log("summary");
                //console.log(summary);
                var summaries = {};

                for (var barcode in summary) {

                    for (var read_type in summary[barcode]) {
                        var lastbinwidth = 0;

                        //Here we make the empty arrays
                        if (summaries[barcode] === undefined) {
                                summaries[barcode] = {};
                            }
                        if (summaries[barcode][read_type] === undefined) {
                                summaries[barcode][read_type] = {
                                    'bin_width': [],
                                    'read_count': [],
                                    'read_length': [],
                                };
                            }

                        var maxbinwidth = 0;

                        for (var tmp_bin_width of Object.keys(summary[barcode][read_type])) {
                            if (parseInt(tmp_bin_width) > maxbinwidth) {
                                maxbinwidth = parseInt(tmp_bin_width);
                            }
                        }

                        //var maxbinwidth = Math.max(Object.keys(summary[barcode][read_type]))

                        //var maxbinwidth = summary[barcode][read_type]
                        for (var i=0;i <= maxbinwidth; i++){
                            summaries[barcode][read_type]['bin_width'].push(parseInt(i * 900 + 900));
                            summaries[barcode][read_type]['read_count'].push(0);
                            summaries[barcode][read_type]['read_length'].push(0);
                        }
                        for (var bin_width in summary[barcode][read_type]) {
                            summaries[barcode][read_type]['read_count'][bin_width]=summary[barcode][read_type][bin_width]['read_count'];
                            summaries[barcode][read_type]['read_length'][bin_width]=summary[barcode][read_type][bin_width]['read_length'];
                            summaries[barcode][read_type]['bin_width'][bin_width]=parseInt(bin_width * 900 + 900);

                        }
                    }
                }


                self.histogramSummary = summaries;
                self.updateHistogramBasedCharts();

            }
        });
    };

    this.updateHistogramBasedCharts = function () {
        self.updateHistogramReadLengthChart();
        self.updateHistogramBasesSequencedReadLengthChart();
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

    this.requestChannelSummaryData = function (id) {
        /*
         * Request channel summary data
         */
        var url = "/api/v1/flowcells/" + id + "/channelsummary_readcount";
        $.get(url, function (data) {
            if (data.length > 0) {
                self.updatePoreChart(self.chartReadsPerPore, data, 'read_count');
            }
        });
        var url = "/api/v1/flowcells/" + id + "/channelsummary_readkb";
        $.get(url, function (data) {
            if (data.length > 0) {
                self.updatePoreChart(self.chartBasesPerPore, data, 'read_length');
            }
        });
    };

    this.updatetext = function () {
        datadump = new Array();
        //console.log('rundetails');
        //console.log(self.rundetails);
        if (self.rundetails != undefined && self.rundetails.length > 0) {
            //console.log(self.rundetails);
            //for (var i = 0; i < self.rundetails.length; i++) {

            for (var i in self.rundetails) {
                // skip loop if the property is from prototype
                if (!self.rundetails.hasOwnProperty(i)) continue;
                //console.log("look at this");
                //console.log(self.rundetails[i]);
                datadump[i] = new Array;
                datadump[i].minIONname = self.rundetails[i]['minION_name'];
                datadump[i].minKNOW_asic_id = self.rundetails[i]['minKNOW_asic_id'];
                datadump[i].minKNOW_version = self.rundetails[i]['minKNOW_version'];
                datadump[i].run_id = self.rundetails[i]['minKNOW_hash_run_id'];
                datadump[i].run_name = self.rundetails[i]['minKNOW_run_name'];
                datadump[i].sample_name = self.rundetails[i]['minKNOW_sample_name'];
                datadump[i].computer_name = self.rundetails[i]['minKNOW_computer'];
                datadump[i].minKNOW_current_script = self.rundetails[i]['minKNOW_current_script'];
                datadump[i].name = self.rundata[i]['name'];

                var starttime = new Date(self.rundetails[i]['minKNOW_start_time']);
                datadump[i].starttime = starttime;
                datadump[i].flowcellid = self.rundetails[i]['minKNOW_flow_cell_id'];
                datadump[i].active = self.rundata[i].active;
                datadump[i].barcodes = self.rundata[i].barcodes;


                //datadump[i].rundetails=self.rundetails[i];
            }
        } else {
            for (var i in self.rundata) {
                datadump[i] = new Array;
                datadump[i].minIONname = self.rundata[i]['minION_name'];
                datadump[i].minKNOW_asic_id = self.rundata[i]['minKNOW_asic_id'];
                datadump[i].minKNOW_version = self.rundata[i]['minKNOW_version'];
                datadump[i].run_id = self.rundata[i]['minKNOW_hash_run_id'];
                datadump[i].run_name = self.rundata[i]['run_name'];
                datadump[i].sample_name = self.rundata[i]['minKNOW_sample_name'];
                datadump[i].computer_name = self.rundata[i]['minKNOW_computer'];
                datadump[i].minKNOW_current_script = self.rundata[i]['minKNOW_current_script'];
                var starttime = new Date(Date.parse(self.rundata[i]['start_time']));
                var lastread = new Date(Date.parse(self.rundata[i]['last_read']));
                datadump[i].starttime = starttime;
                datadump[i].lastread = lastread;
                datadump[i].flowcellid = self.rundata[i]['minKNOW_flow_cell_id'];
                datadump[i].active = self.rundata[i].active;
                datadump[i].barcodes = self.rundata[i].barcodes;

            }
        }
        //console.log("datadump");
        //console.log(datadump);
        var sorteddatadump = datadump.sort(dynamicSort("starttime"));
        self.flowcellstart = sorteddatadump[0].starttime;
        self.write_run_data(sorteddatadump);

    };

    this.requestRunDetails = function (id) {
        var url_RunDetails = '/api/v1/flowcells/' + id + '/rundetails/';
        $.get(url_RunDetails, function (data) {
            //console.log("rundetails");
            //console.log(data);
            self.rundetails = data;
            //self.livedata.minIONname = data[0].minION_name;
            //self.livedata.asicid = data[0].minKNOW_asic_id;
            //self.livedata.scriptid = data[0].minKNOW_current_script;
            if (data[0] != undefined) {
                self.livedata.colours_string = data[0].minKNOW_colours_string;
            }
            //self.livedata.computer=data[0].minKNOW_computer;
            //console.log(data[0]);
            //document.getElementById('ComputerName').innerHTML = data[0].minKNOW_computer;
            //document.getElementById('ScriptID').innerHTML = data[0].minKNOW_current_script;
            self.updatetext();

        })

    };

    this.requestLiveRunStats = function (id) {
        //console.log('lastread ' + this.lastread);
        var url_livestats = '/api/v1/flowcells/' + id + '/runstats/' + this.lastread;
        $.get(url_livestats, function (data) {

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
                self.calculatereadtoeventscaling();

                if (self.needtoupdatecharts == true) {
                    self.updateLiveHistogram(data);
                    self.updateLiveYieldProjection();
                    self.updateLiveCumuYield();
                    self.updatePoreStats();
                    //self.updateTextPredictions();
                }


            }
            //console.log(self.livedata);
        })
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
    };

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
    };

    this.updateLiveCumuYield = function () {
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
    this.converttobases = function (data, seqspeed) {
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
    };

    this.requestMessages = function () {
        //
        var url_sincemessages = self.rundetails[0]["minION"] + 'messagessince/' + self.rundetails[0]["minKNOW_start_time"] + '/' + self.lasttime.toISOString() + "/";
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
    };

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


    this.calculatereadtoeventscaling = function () {
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
    };

    this.requestReference = function (id) {
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
    };
    this.requestTasks = function (id) {
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
    };

    this.liveUpdateTasks = function (id) {
        var url = "/api/v1/flowcells/" + id + "/tasks/";
        $.get(url, function (data) {
            var tasks = [];
            for (var i = 0; i < data.length; i++) {
                tasks.push(data[i]);
            }
            self.tasks = tasks;
            for (var i = 0; i < self.tasks.length; i++) {
                if (self.tasks[i].hasOwnProperty("job_details")) {
                    message = 'Reads processed:' + self.tasks[i]["job_details"]["read_count"] + "/" + self.summary["All reads"]["Template"]["read_count"]["data"][0];
                    percentage = Math.round((self.tasks[i]["job_details"]["read_count"] / self.summary["All reads"]["Template"]["read_count"]["data"][0] * 100) * 100) / 100;
                    message2 = percentage + '% of uploaded reads are processed';
                    //console.log(message);
                    $('#' + self.tasks[i]["name"] + '-message').text(message);
                    $('#' + self.tasks[i]["name"] + '-percentage').text(percentage);
                    $('div#' + self.tasks[i]["name"] + '-percentage').width(percentage + '%');
                    $('#' + self.tasks[i]["name"] + '-message2').text(message2);
                }
            }
        })
    };

    this.updateTasks = function (id) {
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

    };

    this.requestKraken = function (id) {
        var parsedkraken = '/api/v1/flowcells/' + id + '/krakenparse/';
        $.get(parsedkraken, function (data) {
            //console.log(data);
            krakendata = [];
            var list = $("#krakenSelectBarcode");
            var list2 = $("#krakenSelectRead");
            for (var i = 0; i < data.length; i++) {
                var barcodename = data[i]["barcode_name"];
                var readtype = data[i]["type_name"]
                if (!(readtype in krakendata)) {
                    krakendata[readtype] = [];
                    if ($("#krakenSelectRead option[value='" + data[i]["type_name"] + "']").val() === undefined) {
                        list2.append(new Option(data[i]["type_name"], data[i]["type_name"]))
                    }
                }
                if (!(barcodename in krakendata[readtype])) {
                    krakendata[readtype][barcodename] = [];
                    if ($("#krakenSelectBarcode option[value='" + data[i]["barcode_name"] + "']").val() === undefined) {
                        list.append(new Option(data[i]["barcode_name"], data[i]["barcode_name"]));
                    }

                }
                //if (data[i]['percentage'] >= 0.05 && data[i]["parent"] != "Input") {
                //if (data[i]['indentation']>8){
                krakendata[readtype][barcodename].push([data[i]["parent"], data[i]["sci_name"], data[i]['percentage'], data[i]['indentation']])
                //}
            }
            //console.log(krakendata);
            list.change(function () {
                var selectedType = list2.find(":selected").text();
                var selectedBarcode = list.find(":selected").text();
                var minimum = $("#lowerboundselect").val();
                //console.log(minimum);
                var slimmedkraken = [];
                for (var i = 0; i < krakendata[selectedType][selectedBarcode].length; i++) {
                    if (krakendata[selectedType][selectedBarcode][i][3] >= minimum) {
                        slimmedkraken.push([krakendata[selectedType][selectedBarcode][i][0], krakendata[selectedType][selectedBarcode][i][1], krakendata[selectedType][selectedBarcode][i][2]]);
                    }
                }
                ;

                //console.log(selectedBarcode);
                var chart = Highcharts.chart('kraken-sankey', {
                    title: {
                        text: 'Kraken ' + selectedType + ' ' + selectedBarcode
                    },

                    series: [{
                        keys: ['from', 'to', 'weight'],
                        data: slimmedkraken,
                        //data: krakendata[selectedType][selectedBarcode],
                        type: 'sankey',
                        curveFactor: 0,
                        name: 'Kraken Output'
                    }]
                });
            });
        })
    };

    this.requestMappedChromosomes = requestMappedChromosomes;

    this.updateCoverageBasedCharts = function (chart, field) {
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
    };


    this.requestPafData = function (id) {

        self.requestMappedChromosomes(id);

        var pafurl = '/api/v1/flowcells/' + id + '/pafsummary/';

        $.get(pafurl, function (data) {
            //console.log("PAF DATA");
            //console.log(data);
            if (data.length < 1) {
            } else {
                summarycoverage = {};
                for (var i = 0; i < data.length; i++) {
                    //console.log(data[i]);
                    if (summarycoverage[data[i].barcode_group_name] === undefined) {
                        summarycoverage[data[i].barcode_group_name] = {};
                    }
                    if (summarycoverage[data[i].barcode_group_name][data[i].read_type_name] === undefined) {
                        summarycoverage[data[i].barcode_group_name][data[i].read_type_name] = {};
                    }
                    if (summarycoverage[data[i].barcode_group_name][data[i].read_type_name][data[i].chrom_name] === undefined) {
                        summarycoverage[data[i].barcode_group_name][data[i].read_type_name][data[i].chrom_name] = {};
                    }
                    summarycoverage[data[i].barcode_group_name][data[i].read_type_name][data[i].chrom_name]["coverage"] = {
                        "name": "coverage",
                        "data": [data[i].chrom_cover],
                        "animation": false
                    };
                    summarycoverage[data[i].barcode_name][data[i].read_type_name][data[i].chrom_name]["ave_read_len"] = {
                        "name": "Average Read Length",
                        "data": [data[i].avg_read_len],
                        "animation": false
                    };
                }
                self.summarycoverage = summarycoverage;
                self.updateCoverageBasedCharts(self.chart_per_chrom_cov, "coverage");
                self.updateCoverageBasedCharts(self.chart_per_chrom_avg, "ave_read_len");
            }

        })
    }


    this.requestData = function () {

        var url_run = '/api/v1/flowcells/' + self.id;

        $.get(url_run, function (data) {
            //console.log(data);

            var barcodes = new Set();
            for (var j = 0; j < data.length; j++) {
                for (var i = 0; i < data[j].barcodes.length; i++) {
                    barcodes.add(data[j].barcodes[i].name)
                }
            }
            self.barcodes = Array.from(barcodes).sort();
            //console.log(self.barcodes);
            //console.log("This is the data we are logging.");
            //console.log(data);
            self.rundata = data;
            //console.log("rundata");
            //console.log(self.rundata);
            self.updateBarcodeNavTab();
            //self.requestSummaryByMinuteData(self.id);
            self.requestSummaryData(self.id);
            self.requestHistogramData(self.id);
            self.requestChannelSummaryData(self.id);
            self.requestRunDetails(self.id);
            self.requestLiveRunStats(self.id);
            self.liveUpdateTasks(self.id);
            self.requestKraken(self.id);
            self.requestPafData(self.id);
            self.requestCumuBases(self.id);
            self.requestQualTime(self.id);
            self.requestLengthTime(self.id);
            self.requestSpeed(self.id);
            //self.updatetext();


            /*

             //console.log(self.rundata);

             self.requestPafTransData(self.id);

             console.log("seriously - Im just trying to parse kraken");


             */

        });


        console.log(self);

    };

}
