
function afterSetExtremes(e) {

    console.log('---> inside after extremes...');

    var chart = self.chartChromosomeCoverage;

    var select = document.getElementById('chromosome-id-select');
    var selected_index = select.selectedIndex;
    var selected_option = select[selected_index];

    if (selected_option === undefined) {
        console.log('No mapped chromosome selected.');
        return;
    }

    var value_combination = selected_option.value.split('_');
    var run_id = value_combination[0];
    var barcode_id = value_combination[1];
    var read_type_id = value_combination[2];
    var reference_id = value_combination[3];
    var chromosome_id = value_combination[4];


    var url = '/api/v1/runs/' + run_id + '/pafcover/' + barcode_id + '/' + read_type_id + '/' + chromosome_id + '/' + Math.round(e.min) + '/' + Math.round(e.max) + '/';

    //chart.showLoading('Loading data from server...');
    // $.getJSON('https://www.highcharts.com/samples/data/from-sql.php?start=' + Math.round(e.min) +
    //         '&end=' + Math.round(e.max) + '&callback=?', function (data) {
    $.getJSON(url, function (data) {

        var orderedDataSimplified = data['data_simplified'].sort(function (a, b) {
            return a[0] - b[0];
        });

        chart.series[0].setData(orderedDataSimplified);
        chart.hideLoading();
    });
}

function makeStepLineChart(divName, chartTitle, yAxisTitle) {

    var chart = Highcharts.stockChart(divName, {

        chart: {
            zoomType: 'x',
            // panning: true,
            // panKey: 'shift',
            // events: {
            //     selection: function(event) {
            //         var text, label;
            //         if (event.xAxis) {
            //             text = 'min: ' + Highcharts.numberFormat(event.xAxis[0].min, 2) + ', max: ' + Highcharts.numberFormat(event.xAxis[0].max, 2);
            //             self.updateStepLineChart(this, Math.round(event.xAxis[0].min, 0), Math.round(event.xAxis[0].max, 0));
            //         } else {
            //             text = 'Selection reset';
            //             self.updateStepLineChart(this, 0, 0);
            //         }
            //         label = this.renderer.label(text, 100, 120)
            //             .attr({
            //                 fill: Highcharts.getOptions().colors[0],
            //                 padding: 10,
            //                 r: 5,
            //                 zIndex: 8
            //             })
            //             .css({
            //                 color: '#FFFFFF'
            //             })
            //             .add();
            //
            //         setTimeout(function () {
            //             label.fadeOut();
            //         }, 1000);
            //     }
            // }
        },

        navigator: {
            adaptToUpdatedDate: false,
            // series: {
            //     data: data
            // }
        },

        scrollbar: {
            liveRedraw: false
        },

        xAxis: {
            events: {
                afterSetExtremes: self.afterSetExtremes
            }
        },
        // rangeSelector: {
        //     selected: 1
        // },

        title: {
            text: chartTitle
        },
    });

    return chart;

}

function updateStepLineChart (chart, start = 0, end = 0) {

    var select = document.getElementById('chromosome-id-select');
    var selected_index = select.selectedIndex;
    var selected_option = select[selected_index];

    if (selected_option === undefined) {
        console.log('No mapped chromosome selected.');
        return;
    }

    var value_combination = selected_option.value.split('_');
    var run_id = value_combination[0];
    var barcode_id = value_combination[1];
    var read_type_id = value_combination[2];
    var reference_id = value_combination[3];
    var chromosome_id = value_combination[4];

    //var run_id = self.id;
    //var barcode_id = self.selectedBarcode;
    var selected_option_value = selected_option.value;

    var url = '/api/v1/runs/' + run_id + '/pafcover/' + barcode_id + '/' + read_type_id + '/' + chromosome_id + '/' + start + '/' + end + '/';

    $.getJSON(url, function(data) {

        // Remove previous series
        while (chart.series.length > 0) {
            chart.series[0].remove();
        }

        // var orderedDataOriginal = data['data_original'].sort(function (a, b) {
        //     return a[0] - b[0];
        // });

        var orderedDataSimplified = data['data_simplified'].sort(function (a, b) {
            return a[0] - b[0];
        });

        // console.log(orderedData);

        // chart.addSeries({
        //     // name: barcode + " - " + typeName,
        //     name: 'Original',
        //     data: orderedDataOriginal,
        //     step: true,
        //  });

        chart.addSeries({
            // name: barcode + " - " + typeName,
            name: 'Simplified',
            data: orderedDataSimplified,
            step: true,
        });

        console.log('---> orderedDataSimplified...');
        console.log(orderedDataSimplified);

        if (orderedDataSimplified !== undefined) {
            var start = orderedDataSimplified[0][0];
            var end = orderedDataSimplified[orderedDataSimplified.length][0];
            chart.xAxis[0].setExtremes(start, end, true);
        }
    });

}

function requestMappedChromosomes (run_id) {
    /*
     * Request the chromosomes that have reads mapped to using minimap2
     * and update the select box on tab Mapping
     *
     */

    var url = '/api/v1/runs/' + run_id + '/references';

    $.getJSON(url, function(data) {

        var select = document.getElementById('chromosome-id-select');
        var selected_index = select.selectedIndex;
        var selected_option = select[selected_index];

        while (select.length > 0) {
            select.remove(0);
        }

        for (var i = 0; i < data.length; i++) {
            var option = document.createElement('option');

            var value_combination = data[i]['run_id'] + '_' + data[i]['barcode_id'] + '_' + data[i]['read_type_id'] + '_' + data[i]['reference_id'] + '_' + data[i]['chromosome_id'];

            option.text = data[i]['barcode_name'] + ' - ' + data[i]['read_type_name'] + ' - ' + data[i]['reference_name'] + ' - ' + data[i]['chromosome_name'];
            option.value = value_combination;

            if (selected_option === undefined) {

            } else {

                if (value_combination == selected_option.value) {
                    option.selected = 'selected';
                }

            }

            select.add(option);
        }

    });

}
