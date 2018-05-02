function requestHistogramData(id) {
    /*
     * Request histogram data
     */

    var url = "/api/v1/flowcells/" + id + "/histogramsummary/";

    $.get(url, function (data) {

        // console.log(data);

        if (data.length > 0) {

            var selected_barcode = get_selected_barcode();

            var data_barcode = []

            data.forEach(function(row) {

                if (row[0] == selected_barcode) {

                    data_barcode.push(row);
                }
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

            chart.update({
                chart: {
                    type: 'column'
                },
                xAxis: {
                    type: 'category',
                    categories: data_barcode.map(x => x[3])
                }
            });

            chart.addSeries({
                name: data_barcode[0][0] + " - " + data_barcode[0][1] + " - " + data_barcode[0][2],
                data: data_barcode.map(x => x[5])
            });


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
                chart: {
                    type: 'column'
                },
                xAxis: {
                    type: 'category',
                    categories: data_barcode.map(x => x[3])
                }
            });

            chart.addSeries({
                name: data_barcode[0][0] + " - " + data_barcode[0][1] + " - " + data_barcode[0][2],
                data: data_barcode.map(x => x[4])
            });
        }
    });
}
