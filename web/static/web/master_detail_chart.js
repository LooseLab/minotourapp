
function master_detail_after_selection(event) {

    event.preventDefault();

    // console.log(event);
    console.log(event.target.renderTo.id)

    var min = Math.trunc(event.xAxis[0].min);
    var max = Math.trunc(event.xAxis[0].max);

    console.log(event);

    var url = master_detail_create_url() + '?min=' + min + '&max=' + max;

    master_detail_load_chart_data(url);

    master_chart.xAxis[0].removePlotBand('mask-before');
    master_chart.xAxis[0].addPlotBand({
        id: 'mask-before',
        from: min,
        to: max,
        color: 'rgba(0, 0, 0, 0.2)'
    });

}

function master_detail_create_url() {

    return "http://localhost:8000/api/v1/flowcells/1/pafcover/2/1/1/";
}

function master_detail_load_chart_data(url, div_main_name) {

    console.log('before getjson');

    console.log(url);

    $.getJSON(url, function (data) {
        console.log(data);

        console.log('length of master: ' + master_chart.series[0].data.length);
        //console.log('length of detail: ' + detail_chart.series[0].data.length);

        console.log('inside geojson');
        console.log(data);

        if (master_chart.series[0].data.length === 0) {
            master_chart.series[0].setData(data);
        }

        min_extreme = data[0][0];
        max_extreme = data[19][0];
        // min_extreme = data[0][0] + (data[19][0] - data[0][0])/2;
        // max_extreme = data[19][0] - (data[19][0] - data[0][0])/2;

        if (detail_chart) {
            detail_chart.destroy();
        }

        detail_chart = create_detail_chart(div_main_name);

        //detail_chart.xAxis[0].setExtremes(min_extreme, max_extreme);
        detail_chart.series[0].setData(data);
        //detail_chart.series[0].pointStart = data[0][0]

        //detail_chart.redraw();
        return data;

    });

    console.log('after getjson');

}

function create_detail_chart(div_main_name) {

    var div_detail_name = 'teste_div_detail';
    //var div_detail_name = div_main_name + "_detail";
    console.log("detail div name: " + div_detail_name);

    var detail_container = document.getElementById(div_detail_name);

    var chart = Highcharts.chart({
        chart: {
            renderTo: detail_container,
            reflow: false,
            marginLeft: 50,
            marginRight: 20,
            style: {
                position: 'absolute'
            },
            step: true,
            zoomType: 'x',
            events: {
                selection: master_detail_after_selection
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
            // events: {
            //     afterSetExtremes: afterSetExtremes
            // }
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

function create_master_chart(master_div) {

    var chart = new Highcharts.chart({
        chart: {
            renderTo: master_div,
            reflow: false,
            //borderWidth: 0,
            backgroundColor: null,
            marginLeft: 50,
            marginRight: 20,
            zoomType: 'x',
            step: true,
            events: {
                selection: master_detail_after_selection
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

function create_master_detail_chart(div_name) {

    var div_main_name = div_name;
    var div_master_name = div_main_name + "_master";
    // var div_detail_name = div_main_name + "_detail";

    // var detail_container = document.getElementById(div_detail_name);
    var master_container = document.getElementById(div_master_name);
    var main_container = document.getElementById(div_main_name);

    master_chart = create_master_chart(master_container);
    //detail_chart = create_detail_chart(detail_container);

    var url = master_detail_create_url();

    master_detail_load_chart_data(url, div_main_name);

    console.log('done');
}
