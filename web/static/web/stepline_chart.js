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

    $.getJSON('http://localhost:8100/api/v1/runs/6/pafcover/29/1/1/', function(data) {

        var selectedBarcode = self.selectedBarcode;

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

        // for (var barcode of Object.keys(self.summaryByMinute2)) {
        //
        //     if (barcode === 'All reads' || barcode === self.selectedBarcode) {
        //         for (var typeName of Object.keys(self.summaryByMinute2[barcode])) {
        //
        //             chart.addSeries({
        //                 name: barcode + " - " + typeName,
        //                 data: self.summaryByMinute2[barcode][typeName]["sequencingRate"]
        //             });
        //
        //         }
        //     }
        // }

    });

}
