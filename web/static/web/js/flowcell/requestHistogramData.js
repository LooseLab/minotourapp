function requestHistogramData(id) {
    /*
     * Request histogram data
     */

    var url = "/api/v1/flowcells/" + id + "/histogramsummary/";

    $.get(url, function (dataObj) {

        var indexes = dataObj['indexes'];
        var data = dataObj['data'];
        var categories = dataObj['categories'];

        if (data.length > 0) {

            var selected_barcode = get_selected_barcode();

            var data_barcode = []

            data.forEach(function(row) {

                if (row[0] == selected_barcode) {

                    data_barcode.push(row);
                }
            });

            data_barcode.sort(function(a, b) {
                return a[2] - b[2];
            });

            //console.log(data_barcode[-1]);

            indexes.sort(function(a, b) {
                return a - b;
            });

            /*
             * chart
             *
             */
            if (!self.chartHistogramBasesSequencedByReadLength) {
                self.chartHistogramBasesSequencedByReadLength = self.makeChart2(
                    "histogram-bases-sequenced-by-read-length",
                    "Histogram of Bases Sequenced by Read Length".toUpperCase(),
                    "Number of bases".toUpperCase()
                );
            }

            var chart = self.chartHistogramBasesSequencedByReadLength;

            if (chart.series) {
                while (chart.series.length > 0) {
                    chart.series[0].remove();
                }
            }

            //console.log(data_barcode.map(x => x[2]));

            chart.update({
                plotOptions: {
                    column: {
                        stacking: 'normal'
                    }
                },
                chart: {
                    type: 'column'
                },
                xAxis: {
                    type: 'category',
                    categories: data_barcode.map(x => x[2])
                }
            });

            for (var i = 0; i < indexes.length; i++) {

                if (!indexes[i].startsWith("No barcode")) {

                    chart.addSeries({
                        name: indexes[i],
                        data: data_barcode.filter(x => x[1] == indexes[i]).map(x => x[4])
                    });
                }
            }


            /*
             * chart
             *
             */
            if (!self.chartHistogramReadLength) {
                self.chartHistogramReadLength = self.makeChart2(
                    "histogram-read-lengths",
                    "Histogram of Read Lengths".toUpperCase(),
                    "Number of reads".toUpperCase()
                );
            }

            var chart = self.chartHistogramReadLength;

            if (chart.series) {
                while (chart.series.length > 0) {
                    chart.series[0].remove();
                }
            }

            chart.update({
                plotOptions: {
                    column: {
                        stacking: 'normal'
                    }
                },
                chart: {
                    type: 'column'
                },
                xAxis: {
                    type: 'category',
                    categories: data_barcode.map(x => x[2])
                }
            });

            for (var i = 0; i < indexes.length; i++) {

                if (!indexes[i].startsWith("No barcode")) {

                    chart.addSeries({
                        name: indexes[i],
                        data: data_barcode.filter(x => x[1] == indexes[i]).map(x => x[3])
                    });
                }
            }
        }
    });
}
