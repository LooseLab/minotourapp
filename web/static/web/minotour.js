Highcharts.setOptions({
    plotOptions: {
        series: {
            animation: false
        }
    }
});

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


function chartLoadUpdateData(data) {
    var chart_data = data;

    // Remove previous series
    while (this.chart.series.length > 0) {
        this.chart.series[0].remove();
    }

    this.chart.colorCounter = 0;
    this.chart.symbolCounter = 0;

    var summaries = [];
    if (this.mybarcode == 'All reads') {

        for (var i = 0; i < chart_data.length; i++) {
            var average_read_length = chart_data[i].total_length / chart_data[i].read_count;
            var sample_time = new Date(chart_data[i].sample_time);

            var point = {
                x: sample_time,
                y: average_read_length
            }

            summaries.push(point);
        }
    }

    summaries.sort(function (a, b) {
        return a.x - b.x;
    })

    this.chart.addSeries({name: '', data: summaries});
    
}

function makeChart(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: 'column',
            animation: Highcharts.svg, // don't animate in old IE
            marginRight: 10,
        },
        title: {
            text: 'Live random data'
        },
        yAxis: {
            title: {
                text: 'Value'
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
            animation: false
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

function updateChart (chart, data, func) {

}

function MinotourApp() {
    this.chart_reads_called = null;
    this.chart_yield = null;
    this.chart_average_read_length = null;
    this.chart_maximum_read_length = null;
    this.average_read_lengths_overtime = null;

    this.barcodes = null;

    this.summaryByMinute = null;
    this.id = null;
    this.selectedBarcode = null;

    this.makeChart = makeChart;
    this.makeChart2 = makeChart2;

    var self = this;

    this.init = function() {
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

        self.id = document.getElementById('run-id').innerText;
        self.selectedBarcode = 'All Reads';

        this.requestData();

        setInterval(function(){
            this.requestData();
        }.bind(this), 60000);
    }; 

    this.updateChartsBasedOnBarcode = function (event) {
        console.log(event);
        console.log(event.target.innerText);

        self.requestData();
    }

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
            li.appendChild(a);
            ul.appendChild(li);
        }
    }
    
    this.updateSummaryByMinute = function (data) {

    }

    this.requestSummaryByMinuteData = function (id) {
        /*
         * Request summary by minute data
         */
        var url = '/api/v1/runs/' + id + '/summarybyminute';

        $.get(url, function(data) {
            console.log('--- data length: ' + data.length);
            console.log('selected barcador: ' + self.selectedBarcode);

            if (data.length > 0) {
                var orderedData = data.sort(function(a, b){
                    return new Date(a.sample_time) - new Date(b.sample_time);
                });

                /*
                if (self.selectedBarcode == 'All Reads') {
                    self.filteredData = orderedData;
                } else {
                    var filteredData = orderedData.filter(function(el){
                        console.log(el.barcode);
                        console.log(self.selectedBarcode);
                        return el.typename == self.selectedBarcode;
                    });
                }
                */

                self.summaryByMinute = orderedData;

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
            }
        });
    }

    this.requestSummaryData = function (id) {
        /*
         * Request summary by barcode data
         */
        var url = '/api/v1/runs/' + id + '/summarybarcode';

        $.get(url, function(data) {
            console.log('--- data length: ' + data.length);
            console.log('selected barcador: ' + self.selectedBarcode);

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
                            'data': [ item.read_count ], 
                            'animation': false
                        };

                        summaries[item.barcode][item.typename]['yield'] = {
                            'name': item.typename, 
                            'data': [ item.total_length ], 
                            'animation': false
                        };
                        
                        summaries[item.barcode][item.typename]['average_read_length'] = {
                            'name': item.typename,
                            'data': [ item.total_length / item.read_count ],
                            'animation': false
                        };
                        
                        summaries[item.barcode][item.typename]['max_length'] = {
                            'name': item.typename,
                            'data': [ item.max_length ],
                            'animation': false
                        };
                    }
                }

                console.log(summaries);

                while (self.chart_reads_called.series.length  > 0) {
                    self.chart_reads_called.series[0].remove();
                }
                self.chart_reads_called.addSeries(summaries[item.barcode][item.typename]['read_count']);

                while (self.chart_yield.series.length > 0) {
                    self.chart_yield.series[0].remove();
                }
                self.chart_yield.addSeries(summaries[item.barcode][item.typename]['yield']);

                while (self.chart_average_read_length.series.length > 0) {
                    self.chart_average_read_length.series[0].remove();
                }
                self.chart_average_read_length.addSeries(summaries[item.barcode][item.typename]['average_read_length']);

                while (self.chart_maximum_read_length.series.length > 0) {
                    self.chart_maximum_read_length.series[0].remove();
                }
                self.chart_maximum_read_length.addSeries(summaries[item.barcode][item.typename]['max_length']);

            }
        });
    }

    this.requestData = function () {
        console.log('inside request data');

        var url_run = '/api/v1/runs/' + self.id;

        $.get(url_run, function (data) {
            self.barcodes = data.barcodes.sort();
            self.updateBarcodeNavTab();
        });

        self.requestSummaryByMinuteData(self.id);
        self.requestSummaryData(self.id);

        //this.counterDown();
    }
}
