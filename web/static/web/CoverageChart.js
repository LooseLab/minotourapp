function CoverageChart() {

    this.master_chart = null;
    this.detail_chart = null;
}

CoverageChart.prototype.afterSelection = function (event) {

    event.preventDefault();

    var min = Math.trunc(event.xAxis[0].min);
    var max = Math.trunc(event.xAxis[0].max);

    var url = this.create_url() + '?min=' + min + '&max=' + max;

    this.load_chart_data(url);

    this.master_chart.xAxis[0].removePlotBand('mask-before');
    this.master_chart.xAxis[0].addPlotBand({
        id: 'mask-before',
        from: min,
        to: max,
        color: 'rgba(0, 0, 0, 0.2)'
    });

}

CoverageChart.prototype.create_url = function () {

    var selected_index = this.select_container.selectedIndex;
    var selected_option = this.select_container[selected_index];

    if (parseInt(selected_option.value) < 0) {
        console.log('No mapped chromosome selected.');
        return;
    }

    console.log('Selected option: ' + selected_option.value);

    // var value_combination = data[i]['task_id'] + '_' + data[i]['barcode_name'] + '_' + data[i]['read_type_id'] + '_' + '_' + data[i]['chromosome_id'];

    var value_combination = selected_option.value.split('_');
    var task_id = value_combination[0];
    // var run_id = value_combination[0];
    var barcode_name = value_combination[1];
    // var barcode_id = value_combination[1];
    var read_type_id = value_combination[2];
    // var reference_id = value_combination[3];
    var chromosome_id = value_combination[3];
    // var chromosome_id = value_combination[4];

    //return "http://localhost:8000/api/v1/flowcells/1/pafcover/2/1/1/";
    var url = "/api/v1/flowcells/" + task_id + "/pafcover/" + barcode_name + "/" + read_type_id + "/" + chromosome_id + "/";
    // var url = "/api/v1/flowcells/" + run_id + "/pafcover/" + barcode_id + "/" + read_type_id + "/" + chromosome_id + "/";

    return url;
}

CoverageChart.prototype.load_chart_data = function (url, div_main_name) {

    $.getJSON(url, (function (data) {

        console.log('>>> got data from server');
        console.log('master_chart series data length is: ' + this.master_chart.series[0].data.length);

        if (this.master_chart.series[0].data.length === 0) {
            console.log('setting master_chart data');
            this.master_chart.series[0].setData(data);
        }

        min_extreme = data[0][0];
        max_extreme = data[19][0];
        // min_extreme = data[0][0] + (data[19][0] - data[0][0])/2;
        // max_extreme = data[19][0] - (data[19][0] - data[0][0])/2;

        if (this.detail_chart) {
            this.detail_chart.destroy();
        }

        this.detail_chart = this.create_detail_chart(div_main_name);

        //detail_chart.xAxis[0].setExtremes(min_extreme, max_extreme);
        this.detail_chart.series[0].setData(data);
        //detail_chart.series[0].pointStart = data[0][0]

        //detail_chart.redraw();
        return data;

    }).bind(this));
}

CoverageChart.prototype.create_detail_chart = function () {

    var chart = Highcharts.chart({
        chart: {
            renderTo: this.detail_container,
            reflow: false,
            marginLeft: 50,
            marginRight: 20,
            style: {
                position: 'absolute'
            },
            step: true,
            zoomType: 'x',
            events: {
                selection: this.afterSelection.bind(this)
            }
        },
        credits: {
            enabled: false
        },
        title: {
            text: null
        },
        xAxis: {
            type: 'line',
        },
        yAxis: {
            title: {
                text: null
            },
            maxZoom: 0.1,
            min: 0
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            series: {
                marker: {
                    enabled: false,
                    states: {
                        hover: {
                            enabled: true,
                            radius: 3
                        }
                    }
                }
            }
        },
        series: [{
            name: null,
            //pointStart: detailStart,
            data: [],
            step: true
        }],

        exporting: {
            enabled: false
        }

    });

    return chart;

}

CoverageChart.prototype.create_master_chart = function () {

    var chart = new Highcharts.chart({
        chart: {
            renderTo: this.master_container,
            reflow: true,
            //borderWidth: 0,
            backgroundColor: null,
            marginLeft: 50,
            marginRight: 20,
            zoomType: 'x',
            step: true,
            events: {
                selection: this.afterSelection.bind(this)
            }
        },
        title: {
            text: null
        },
        xAxis: {
            type: 'line',
            showLastTickLabel: true,
            plotBands: [{
                id: 'mask-before',
                //from: data[0][0],
                //to: data[data.length - 1][0],
                color: 'rgba(0, 0, 0, 0.2)'
            }],
            title: {
                text: null
            }
        },
        yAxis: {
            title: {
                text: null
            },
            visible: false
        },
        legend: {
            enabled: false
        },
        credits: {
            enabled: false
        },
        plotOptions: {
            series: {
                fillColor: {
                    // linearGradient: [0, 0, 0, 70],
                    stops: [
                        [0, Highcharts.getOptions().colors[0]],
                        [1, 'rgba(255,255,255,0)']
                    ]
                },
                lineWidth: 1,
                marker: {
                    enabled: false
                },
                shadow: false,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                enableMouseTracking: false
            }
        },

        series: [{
            type: 'area',
            name: null,
            //pointStart: data[0][0],
            data: [],
            step: true
        }],

        exporting: {
            enabled: false
        }

    });

    return chart;

}

CoverageChart.prototype.create_master_detail_chart = function (div_name, select_id) {

    var div_main_name = div_name;

    this.div_master_name = div_main_name + "_master";
    this.div_detail_name = div_main_name + "_detail";

    this.master_container = document.getElementById(this.div_master_name);
    this.detail_container = document.getElementById(this.div_detail_name);

    this.master_chart = this.create_master_chart();

    this.select_container = document.getElementById(select_id);

    this.select_container.addEventListener('change', (function(event){

        var url = this.create_url();

        if (url) {

            this.load_chart_data(url, div_main_name);
        }

        //$('#chromosome-id-select').on('change', function () {
        //    self.updateStepLineChart(self.chartChromosomeCoverage, 0, 0);
        //});
    }).bind(this));
}
