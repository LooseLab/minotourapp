function makeStepLineChart(divName, chartTitle, yAxisTitle) {

    var chart = Highcharts.chart(divName, {

        chart: {
            zoomType: 'x'
        },

        rangeSelector: {
            selected: 1
        },

        title: {
            text: chartTitle
        },
    });

    return chart;

}

function updateStepLineChart (chart) {

    var select = document.getElementById('chromosome-id-select');
    var selected_index = select.selectedIndex;
    var selected_option = select[selected_index];

    var value_combination = selected_option.value.split('_');
    var run_id = value_combination[0];
    var barcode_id = value_combination[1];
    var read_type_id = value_combination[2];
    var reference_id = value_combination[3];
    var chromosome_id = value_combination[4];

    if (selected_option === undefined) {
        console.log('No mapped chromosome selected.');
        return;
    }

    //var run_id = self.id;
    //var barcode_id = self.selectedBarcode;
    var selected_option_value = selected_option.value;

    var url = '/api/v1/runs/' + run_id + '/pafcover/' + barcode_id + '/' + read_type_id + '/' + chromosome_id + '/';

    $.getJSON(url, function(data) {

        // Remove previous series
        while (chart.series.length > 0) {
            chart.series[0].remove();
        }

        chart.addSeries({
            // name: barcode + " - " + typeName,
            name: 'Test',
            data: data,
            step: true,
         });

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
