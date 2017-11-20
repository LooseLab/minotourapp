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

function updateStepLineChart (chart, run_id, barcode_id) {

    var select = document.getElementById('chromosome-id-select');
    var selected_index = select.selectedIndex;
    var selected_option = select[selected_index];

    if (selected_option === undefined) {
        console.log('No mapped chromosome selected.');
        return;
    }

    //var run_id = self.id;
    //var barcode_id = self.selectedBarcode;
    var selected_option_value = selected_option.value;

    var url = '/api/v1/runs/' + run_id + '/pafcover/' + barcode_id + '/1/' + selected_option_value + '/';

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

            option.text = data[i]['barcode__name'] + ' - ' + data[i]['reference_name'] + ' - ' + data[i]['chromosome__name'];
            option.value = data[i]['chromosome_id'];

            if (selected_option === undefined) {

            } else {

                if (data[i][0] == selected_option.value) {
                    option.selected = 'selected';
                    console.log('selected option is ' + selected_option.value);
                }

            }

            select.add(option);
        }

    });

}