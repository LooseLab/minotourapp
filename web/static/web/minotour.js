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


function MinotourApp() {
    this.chart_reads_called = null;
    this.chart_yield = null;
    this.chart_average_read_length = null;
    this.chart_maximum_read_length = null;
    this.average_read_lengths_overtime = null;
    this.chart_chart_cumulative_number_reads_overtime = null;
    this.chartSequencingRate = null;

    this.barcodes = null;

    this.summaryByMinute = null;
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
    }

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
    }

    this.updateCumulativeNumberOfReadsOverTime = function() {
        // chart chart_cumulative_number_reads_overtime
        var summaries = {};

        for (var i = 0; i < self.summaryByMinute.length; i++) {
            if (summaries[self.summaryByMinute[i].typename] === undefined) {
                summaries[self.summaryByMinute[i].typename] = {
                    'lastCumulativeReadCount': 0,
                    'data': []
                };
            }

            var cumulativeReadCount  = self.summaryByMinute[i].read_count + summaries[self.summaryByMinute[i].typename].lastCumulativeReadCount;
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

    }

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
        self.chartSequencingRate.update({
            chart: {
                type: 'column'
            }
        });
    }

    this.updateCumulativeNumberOfReadsOverTime = function() {

        // Remove previous series
        while (self.chart_cumulative_number_reads_overtime.series.length > 0) {
            self.chart_cumulative_number_reads_overtime.series[0].remove();
        }

        var summaries = {};

        for (var i = 0; i < self.summaryByMinute.length; i++) {
            if (summaries[self.summaryByMinute[i].typename] === undefined) {
                summaries[self.summaryByMinute[i].typename] = {
                    'lastCumulativeReadCount': 0,
                    'data': []
                };
            }

            var cumulativeReadCount  = self.summaryByMinute[i].read_count + summaries[self.summaryByMinute[i].typename].lastCumulativeReadCount;
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

    }


    this.updateSummaryByMinuteBasedCharts = function () {
        self.updateCumulativeNumberOfReadsOverTime();
        self.updateSequencingRateChart();
    }

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

                var filteredData = orderedData.filter(function (el) {
                    return el.barcode === self.selectedBarcode;
                });

                self.summaryByMinute = filteredData;

                // Remove previous series
                while (self.average_read_lengths_overtime.series.length > 0) {
                    self.average_read_lengths_overtime.series[0].remove();
                }

                var summaries = {};

                for (var i = 0; i < self.summaryByMinute.length; i++) {
                    if (summaries[self.summaryByMinute[i].typename] === undefined) {
                        summaries[self.summaryByMinute[i].typename] = [];
                    }
                    var average_read_length = self.summaryByMinute[i].total_length / self.summaryByMinute[i].read_count;
                    var sample_time = new Date(self.summaryByMinute[i].sample_time);

                    var point = {
                        x: sample_time,
                        y: average_read_length
                    }

                    summaries[self.summaryByMinute[i].typename].push(point);
                }

                for (var key in summaries) {
                    self.average_read_lengths_overtime.addSeries({name: key, data: summaries[key]});
                }

                // update chart - TODO split ajax call and update of self.summaryByMinute from the redrawing the charts
                self.updateSummaryByMinuteBasedCharts();
            }
        });
    }

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

                var charts = {
                    'read_count': self.chart_reads_called,
                    'yield': self.chart_yield,
                    'average_read_length': self.chart_average_read_length,
                    'max_length': self.chart_maximum_read_length
                };


                for (var prop in charts) {
                    while (charts[prop].series.length > 0) {
                        charts[prop].series[0].remove();
                    }

                    for (key in summaries[self.selectedBarcode]) {
                        charts[prop].addSeries(summaries[self.selectedBarcode][key][prop]);
                    }
                }
            }
        });
    }

    this.requestData = function () {
        var url_run = '/api/v1/runs/' + self.id;

        $.get(url_run, function (data) {
            //console.log(data);
            self.barcodes = data.barcodes.sort();
            self.updateBarcodeNavTab();
        });

        self.requestSummaryByMinuteData(self.id);
        self.requestSummaryData(self.id);
    }
}
