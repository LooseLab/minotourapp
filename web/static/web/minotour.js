Highcharts.setOptions({
    plotOptions: {
        series: {
            animation: false,
            turboThreshold: 0,
        }
    }
});

var NUMBER_SECONDS_IN_A_MINUTE = 60;

var MINOTOUR_VERSION = 0.5;

function check_minotour_version() {
    $.getJSON('http://www.nottingham.ac.uk/~plzloose/minoTourhome/message.php?callback=?', function (result) {

        $.each(result, function (key, value) {
            //checking version info.
            if (key == 'version') {
                if (value == MINOTOUR_VERSION) {
                    $('#newstarget').html("You are running the most recent version of minoTour - version " + value + ".<br>");
                }
                else if (value < MINOTOUR_VERSION) {
                    $('#newstarget').html("You appear to be in the fortunate position of running a future version of the minoTour web application " + value + ". If you have modified the code yourself - great. If not then there might be an issue somewhere!.<br>");
                }
                else if (value > MINOTOUR_VERSION) {
                    $('#newstarget').html("You are running an outdated version of the minoTour web application. The most recent version of minoTour is version " + value + ".<br>" + "Instructions for upgrading will be posted below.<br>");
                }

            } else if (key.substring(0, 7) == 'message') {
                $('#newstarget').append(value + "<br>");
            }
        });
    });
}

function check_user_runs() {
    var url = '/api/v1/runs/';

    $.getJSON(url, function (data) {
        var items = [];
        $.each(data, function (key, val) {
            console.log(val);
            items.push("<li id='" + key + "'><a href='/web/private/runs/" + val.id + "'>" + val.run_name + "</a></li>");
        });

        var text = "<span>You have " + data.length + " minION runs available to view.</span>"
        text += "<ul>";
        text += items.join("");
        text += "</ul>";

        $("#div-user-runs").html(text);
    });
}


function makeChart(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: 'column',
            animation: Highcharts.svg, // don't animate in old IE
            marginRight: 10,
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
                color: '#808080'
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
            type: 'spline',
            marginRight: 10,
            animation: false,
            zoomType: 'x',
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
                text: yAxisTitle
            },
            plotLines: [{
                value: 0,
                width: 1,
                color: '#808080'
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

function MonitorAPP() {
    this.livedata = new Array();
    var self = this;

    this.init = function () {
        console.log('This is MonitorApp Running');
        this.requestData();

        setInterval(function () {
            this.requestData();
        }.bind(this), 10000);
    }

    this.updatecounters = function (live) {
        console.log('update called');
        var thing = document.getElementById('livenum');
        thing.innerHTML = live;
        var thing2 = document.getElementById('livenumruns');
        thing2.innerHTML = 'You have ' + live + ' live runs.';
    }

    this.requestData = function () {
        var url_run = '/api/v1/currentruns/';

        $.get(url_run, function (data) {
            //console.log(data);
            self.livedata = data;
            self.updatecounters(self.livedata.length + '/' + '0');
            //self.barcodes = data.barcodes.sort();
            //self.updateBarcodeNavTab();
        }.bind(this));
        console.log(self.livedata.length);
        //self.updatecounters(self.livedata.length + '/' + '0');
        //self.requestSummaryByMinuteData(self.id);
        //self.requestSummaryData(self.id);
    }
}

function MinotourApp() {
    console.log('This is MinotourApp Running');
    this.chart_reads_called = null;
    this.chart_yield = null;
    this.chart_average_read_length = null;
    this.chart_maximum_read_length = null;
    this.average_read_lengths_overtime = null;
    this.chart_chart_cumulative_number_reads_overtime = null;
    this.chartSequencingRate = null;

    this.barcodes = null;

    this.summaryByMinute = null;
    this.summary = null;

    this.id = null;
    this.selectedBarcode = null;

    this.makeChart = makeChart;
    this.makeChart2 = makeChart2;

    var self = this;

    this.init = function () {
        this.chart_reads_called = this.makeChart(
            'reads-called',
            'reads called'.toUpperCase(),
            'number of reads called'.toUpperCase()
        );

        this.chart_yield = this.makeChart(
            'yield',
            'yield'.toUpperCase(),
            'yield'.toUpperCase()
        );

        this.chart_average_read_length = this.makeChart(
            'average-read-length',
            'average read length'.toUpperCase(),
            'average read length'.toUpperCase()
        );

        this.chart_maximum_read_length = this.makeChart(
            'maximum-read-length',
            'maximum read length'.toUpperCase(),
            'maximum read length'.toUpperCase()
        );

        this.average_read_lengths_overtime = this.makeChart2(
            'average-read-lengths-overtime',
            'average read length over time'.toUpperCase(),
            'average read length'.toUpperCase()
        );

        this.chart_cumulative_number_reads_overtime = this.makeChart2(
            'cumulative-number-reads-overtime',
            'cumulative reads'.toUpperCase(),
            'cumulative reads'.toUpperCase()
        );

        this.chartSequencingRate = this.makeChart2(
            'sequencing-rate',
            'sequencing rate'.toUpperCase(),
            'bases/second'.toUpperCase()
        );


        self.id = document.getElementById('run-id').innerText;
        self.selectedBarcode = 'All reads';

        this.requestData();

        setInterval(function () {
            this.requestData();
        }.bind(this), 10000);
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
        var ul = document.getElementById('nav-tabs-barcodes');

        ul.innerHTML = "";

        var sortedBarcodes = this.barcodes;

        for (var i = 0; i < sortedBarcodes.length; i++) {
            var li = document.createElement('li');
            var a = document.createElement('a');
            a.onclick = this.updateChartsBasedOnBarcode;
            a.href = '#';
            a.text = sortedBarcodes[i];

            if (sortedBarcodes[i] === self.selectedBarcode) {
                li.classList.add('active');
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

        for (var readtype of Object.keys(summaries['All reads'])) {
            data.push(summaries['All reads'][readtype][field]['data'][0]);
        }

        serie = {
            'name': 'All reads',
            'data': data
        };

        series.push(serie);

        // Include specific barcode if selected
        if (self.selectedBarcode !== 'All reads') {

            data = [];

            for (var readtype of Object.keys(summaries[self.selectedBarcode])) {
                data.push(summaries[self.selectedBarcode][readtype][field]['data'][0]);
            }

            serie = {
                'name': self.selectedBarcode,
                'data': data
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

        if (selectedBarcode !== 'All reads') {
            var summaries = {};
            summaries['All reads'] = {};
            summaries[selectedBarcode] = {};

        } else {
            var summaries = {
                'All reads': {}
            };

        }

        for (var i = 0; i < self.summaryByMinute.length; i++) {

            var average_read_length = self.summaryByMinute[i].total_length / self.summaryByMinute[i].read_count;
            var sample_time = new Date(self.summaryByMinute[i].sample_time);

            var point = {
                x: sample_time,
                y: average_read_length
            }

            if (self.summaryByMinute[i].barcode === 'All reads') {
                if (summaries['All reads'][self.summaryByMinute[i].typename] === undefined) {
                    summaries['All reads'][self.summaryByMinute[i].typename] = [];
                }

                summaries['All reads'][self.summaryByMinute[i].typename].push(point);
            }

            if (self.summaryByMinute[i].barcode === selectedBarcode && selectedBarcode !== 'All reads') {
                if (summaries[selectedBarcode][self.summaryByMinute[i].typename] === undefined) {
                    summaries[selectedBarcode][self.summaryByMinute[i].typename] = [];
                }

                summaries[selectedBarcode][self.summaryByMinute[i].typename].push(point);
            }

        }

        for (var barcode in summaries) {
            for (var readtype in summaries[barcode]) {
                chart.addSeries({name: barcode + ' - ' + readtype, data: summaries[barcode][readtype]});
            }
        }
    };

    this.updateCumulativeNumberOfReadsOverTimeChart = function () {
        var chart = self.chart_cumulative_number_reads_overtime;
        var selectedBarcode = self.selectedBarcode;

        while (chart.series.length > 0) {
            chart.series[0].remove();
        }

        if (selectedBarcode !== 'All reads') {
            var summaries = {};
            summaries['All reads'] = {};
            summaries[selectedBarcode] = {};

        } else {
            var summaries = {
                'All reads': {}
            };

        }

        for (var i = 0; i < self.summaryByMinute.length; i++) {

            if (self.summaryByMinute[i].barcode === 'All reads') {
                if (summaries['All reads'][self.summaryByMinute[i].typename] === undefined) {
                    summaries['All reads'][self.summaryByMinute[i].typename] = {
                        'lastCumulativeReadCount': 0,
                        'data': []
                    };
                }

                var cumulativeReadCount = self.summaryByMinute[i].read_count + summaries['All reads'][self.summaryByMinute[i].typename].lastCumulativeReadCount;
                var sample_time = new Date(self.summaryByMinute[i].sample_time);

                var point = {
                    x: sample_time,
                    y: cumulativeReadCount
                }

                summaries['All reads'][self.summaryByMinute[i].typename].lastCumulativeReadCount = cumulativeReadCount;
                summaries['All reads'][self.summaryByMinute[i].typename].data.push(point);
            }

            if (self.summaryByMinute[i].barcode === selectedBarcode && selectedBarcode !== 'All reads') {
                if (summaries[selectedBarcode][self.summaryByMinute[i].typename] === undefined) {
                    summaries[selectedBarcode][self.summaryByMinute[i].typename] = {
                        'lastCumulativeReadCount': 0,
                        'data': []
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

        console.log(summaries);
        for (var barcode in summaries) {
            for (var readtype in summaries[barcode]) {
                chart.addSeries({name: barcode + ' - ' + readtype, data: summaries[barcode][readtype]['data']});
            }
        }
    };

    this.updateCumulativeNumberOfReadsOverTimeChart2 = function () {
        // chart chart_cumulative_number_reads_overtime
        var summaries = {};

        for (var i = 0; i < self.summaryByMinute.length; i++) {
            if (summaries[self.summaryByMinute[i].typename] === undefined) {
                summaries[self.summaryByMinute[i].typename] = {
                    'lastCumulativeReadCount': 0,
                    'data': []
                };
            }

            var cumulativeReadCount = self.summaryByMinute[i].read_count + summaries[self.summaryByMinute[i].typename].lastCumulativeReadCount;
            summaries[self.summaryByMinute[i].typename].lastCumulativeReadCount = cumulativeReadCount;

            var sample_time = new Date(self.summaryByMinute[i].sample_time);

            var point = {
                x: sample_time,
                y: cumulativeReadCount
            }

            summaries[self.summaryByMinute[i].typename]['data'].push(point);
        }

        for (var key in summaries) {
            self.chart_cumulative_number_reads_overtime.addSeries({name: key, data: summaries[key]['data']});
        }

    };

    this.updateSequencingRateChart = function () {

        // Remove previous series
        while (self.chartSequencingRate.series.length > 0) {
            self.chartSequencingRate.series[0].remove();
        }

        var summaries = {};

        for (var i = 0; i < self.summaryByMinute.length; i++) {
            if (summaries[self.summaryByMinute[i].typename] === undefined) {
                summaries[self.summaryByMinute[i].typename] = {
                    'data': []
                };
            }

            var sequencingRate = self.summaryByMinute[i].total_length / NUMBER_SECONDS_IN_A_MINUTE;

            var sample_time = new Date(self.summaryByMinute[i].sample_time);

            var point = {
                x: sample_time,
                y: sequencingRate
            }

            summaries[self.summaryByMinute[i].typename]['data'].push(point);
        }

        for (var key in summaries) {
            self.chartSequencingRate.addSeries({name: key, data: summaries[key]['data']});
        }
    };

    this.updateSummaryByMinuteBasedCharts = function () {
        self.updateCumulativeNumberOfReadsOverTimeChart();
        self.updateSequencingRateChart();
        self.updateAverageReadLengthOverTimeChart()
    };

    this.updateSummaryBasedCharts = function () {
        var charts = {
            'read_count': self.chart_reads_called,
            'yield': self.chart_yield,
            'average_read_length': self.chart_average_read_length,
            'max_length': self.chart_maximum_read_length
        };

        for (var prop in charts) {
            self.updateReadsColumnBasedChart(charts[prop], prop);
        }
    };

    this.requestSummaryByMinuteData = function (id) {
        /*
         * Request summary by minute data
         */
        var url = '/api/v1/runs/' + id + '/summarybarcodebyminute';

        $.get(url, function (data) {

            if (data.length > 0) {
                var orderedData = data.sort(function (a, b) {
                    return new Date(a.sample_time) - new Date(b.sample_time);
                });

                self.summaryByMinute = orderedData;

                // update chart - TODO split ajax call and update of self.summaryByMinute from the redrawing the charts
                self.updateSummaryByMinuteBasedCharts();
            }
        });
    };

    this.requestSummaryData = function (id) {
        /*
         * Request summary by barcode data
         */
        var url = '/api/v1/runs/' + id + '/summarybarcode';

        $.get(url, function (data) {

            if (data.length > 0) {
                var summaries = {};
                for (var i = 0; i < data.length; i++) {
                    var item = data[i];

                    if (summaries[item.barcode] === undefined) {
                        summaries[item.barcode] = {};
                    }

                    if (summaries[item.barcode][item.typename] === undefined) {
                        summaries[item.barcode][item.typename] = {
                            'read_count': null,
                            'yield': null,
                            'average_read_length': null,
                            'max_length': null
                        };

                        summaries[item.barcode][item.typename]['read_count'] = {
                            'name': item.typename,
                            'data': [item.read_count],
                            'animation': false
                        };

                        summaries[item.barcode][item.typename]['yield'] = {
                            'name': item.typename,
                            'data': [item.total_length],
                            'animation': false
                        };

                        summaries[item.barcode][item.typename]['average_read_length'] = {
                            'name': item.typename,
                            'data': [item.total_length / item.read_count],
                            'animation': false
                        };

                        summaries[item.barcode][item.typename]['max_length'] = {
                            'name': item.typename,
                            'data': [item.max_length],
                            'animation': false
                        };
                    }
                }

                self.summary = summaries;

                self.updateSummaryBasedCharts();

            }
        });
    };

    this.requestData = function () {
        var url_run = '/api/v1/runs/' + self.id;

        $.get(url_run, function (data) {
            self.barcodes = data.barcodes.sort();
            self.updateBarcodeNavTab();
        });

        self.requestSummaryByMinuteData(self.id);
        self.requestSummaryData(self.id);
    };
}
