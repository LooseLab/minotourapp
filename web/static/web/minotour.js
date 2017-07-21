import jQuery from 'jquery';
import Highcharts from 'highcharts';
//import Vue from 'vue';
//import VueResource from 'vue-resource';

// https://stackoverflow.com/questions/34338411/how-to-import-jquery-using-es6-syntax
window.$ = window.jQuery = jQuery;
//window.Vue = Vue;

//Vue.use(VueResource);

require('bootstrap');
//require('highcharts');
require('datatables');

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

function chart_cumulative_reads_init() {
    console.log('inside chart_cumulative_reads_init.')
}

function chart_cumulative_reads_update() {
    console.log('inside chart_cumulative_reads_update.')

    if (this.chart_data == null) {
        console.log('---- chart_data is null ----');

    } else {
        console.log('---- chart_data is NOT null ----');

        var chart_data = this.chart_data;

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

}

function chartLoadInitialData(data) {

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

function makeChart(divName) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: 'spline',
            animation: Highcharts.svg, // don't animate in old IE
            marginRight: 10,
            events: {
                load: function () {

                    // set up the updating of the chart each second
                    var series = this.series[0];
                    setInterval(function () {
                        var x = (new Date()).getTime(), // current time
                            y = Math.random();
                        series.addPoint([x, y], true, true);
                    }, 1000);
                }
            }
        },
        title: {
            text: 'Live random data'
        },
        xAxis: {
            type: 'datetime',
            tickPixelInterval: 150
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
        tooltip: {
            formatter: function () {
                return '<b>' + this.series.name + '</b><br/>' +
                    Highcharts.dateFormat('%Y-%m-%d %H:%M:%S', this.x) + '<br/>' +
                    Highcharts.numberFormat(this.y, 2);
            }
        },
        legend: {
            enabled: false
        },
        exporting: {
            enabled: false
        },
        series: [{
            name: 'Random data',
            data: (function () {
                // generate an array of random data
                var data = [],
                    time = (new Date()).getTime(),
                    i;

                for (i = -19; i <= 0; i += 1) {
                    data.push({
                        x: time + i * 1000,
                        y: Math.random()
                    });
                }
                return data;
            }())
        }]
    });

    return chart;
}


function makeChart2(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: 'spline',
            animation: Highcharts.svg, // don't animate in old IE
            marginRight: 10,
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


function test() {
    return 3;
}







function updateChartReadsCalled() {
    if (this.chart_reads_called) {

    }

}

function updateChartsBasedOnBarcode (event) {
    console.log(event);
    console.log(event.target.innerText);

    this.updateChartReadsCalled();

}

function updateBarcodeNavTab () {
    var ul = document.getElementById('nav-tabs-barcodes');
    var sortedBarcodes = this.barcodes();

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

function requestInitialData () {
    console.log('inside request initial data');

    const id = 1;

    var url_run = '/api/v1/runs/' + id;

    $.get(url_run, function (data) {
        console.log(data);
        console.log(data.barcodes);

        this.barcodes = data.barcodes.sort();
        this.updateBarcodeNavTab();
    });

    /*
     .then(function (minion_run) {
     this.selectedBarcode = 'All reads';
     this.minion_run = minion_run;
     this.minion_run.barcodes = this.minion_run.barcodes.sort()
     });

     url_run = '/api/v1/runs/' + this.id + '/summarybyminute';
     //promise_run = this.$http.get(url_run);

     promise_run =
     */
    /*
     promise_run
     .then(function (res) {
     return res.json();
     })
     .then(function (summarybyminute) {
     this.summary_by_minute = summarybyminute;
     console.log('---- summary by minute ----');
     console.log(summarybyminute);
     console.log('---- summary by minute ----');
     });

     //this.counterDown();
     */
}

function MinotourApp() {
    this.chart_reads_called = null;
    this.chart_yield = null;
    this.chart_average_read_length = null;
    this.chart_maximum_read_length = null;
    this.average_read_lengths_overtime = null;

    this.barcodes = null;

    this.summaryByMinute = null;

    this.updateChartsBasedOnBarcode = updateChartsBasedOnBarcode;
    this.makeChart = makeChart;
    this.makeChart2 = makeChart2;

    var self = this;

    this.init = function() {
        this.chart_reads_called = this.makeChart('reads-called');
        this.chart_yield = this.makeChart('yield');
        this.chart_average_read_length = this.makeChart('average-read-length');
        this.chart_maximum_read_length = this.makeChart('maximum-read-length');
        this.average_read_lengths_overtime = this.makeChart2(
            'average-read-lengths-overtime', 
            'average read length over time',
            'average read length'.toUpperCase
        );

        //this.request
        setInterval(function(){
            this.test();
            this.requestData();
        }.bind(this), 5000);
    }; 

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
    
    this.requestInitialData = function () {
        console.log('inside request initial data');

        const id = 1;

        var url_run = '/api/v1/runs/' + id;

        $.get(url_run, function (data) {
            console.log(data);
            console.log(data.barcodes);

            self.barcodes = data.barcodes.sort();
            self.updateBarcodeNavTab();
        });

        $.get('/api/v1/runs/' + id + '/summarybyminute', function(data) {
            self.summaryByMinute = data;
            console.log('---- summary by minute ----');
            console.log(self.summaryByMinute);
            console.log('---- summary by minute ----');
        });

        //this.counterDown();
    }

    this.updateSummaryByMinute = function (data) {

    }

    this.requestData = function () {
        console.log('inside request data');

        const id = 1;

        var url_run = '/api/v1/runs/' + id;

        $.get(url_run, function (data) {
            self.barcodes = data.barcodes.sort();
            self.updateBarcodeNavTab();
        });

        /*
         * Request summary by minute data
         */
        var url = '';
        if (!this.summaryByMinute) {
            console.log('-- summaries by minute equal null.');
            url = '/api/v1/runs/' + id + '/summarybyminute';

        } else {
            console.log('-- summaries by minute NOT equal null.');
            var lastSummary = this.summaryByMinute[this.summaryByMinute.length - 1];
            console.log(lastSummary.sample_time);
            url = '/api/v1/runs/' + id + '/summarybyminute/' + lastSummary.sample_time;
        }

        $.get(url, function(data) {
            console.log('--- data length: ' + data.length);

            if (data.length > 0) {
                var orderedData = data.sort(function(a, b){
                    return new Date(a.sample_time) - new Date(b.sample_time);
                });
                
                console.log(orderedData);

                if (self.summaryByMinute) {
                    orderedData.forEach(function(element) {
                        self.summaryByMinute.push(element);

                        self.average_read_lengths_overtime.series.forEach(function(theSerie){
                            var average_read_length = element.total_length / element.read_count;
                            var sample_time = new Date(element.sample_time);

                            var point = {
                                x: sample_time,
                                y: average_read_length
                            }

                            if (theSerie.name === element.typename) {
                                theSerie.addPoint(point);
                            } else {
                                self.average_read_lengths_overtime.addSeries({name: element.typename, data: [point]});
                            }
                        })
                    }, this);
                } else {
                    self.summaryByMinute = orderedData;

                    if (self.average_read_lengths_overtime.series.length == 0) {
                        console.log(self.average_read_lengths_overtime.series);

                        var summaries = {};
                        //if (this.mybarcode == 'All reads') {
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
                        //}

                        //summaries.sort(function (a, b) {
                        //    return a.x - b.x;
                        //})

                        for (var key in summaries) {
                            self.average_read_lengths_overtime.addSeries({name: key, data: summaries[key]});
                        }
                        


                    }



                }
            }
        });

        //this.counterDown();
    }

    this.test = function() {
        console.log('teste');
    }
}







var app = new MinotourApp();

//app.requestInitialData();
app.init();

