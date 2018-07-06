function requestHistogramData(id) {
    /*
     * Request histogram data
     */

    var selected_barcode = get_selected_barcode();

    var url = "/api/v1/flowcells/" + id + "/histogramsummary/";

    $.get(url, function (dataObj) {

        var result_read_count_sum = dataObj['result_read_count_sum'];
        var result_read_length_sum = dataObj['result_read_length_sum'];
        var categories = dataObj['categories'];

        if (data.length > 0) {

            var selected_barcode = get_selected_barcode();

            var data_barcode = []

            var categories = []

            var category_labels = []

            data.forEach(function(row) {
                //console.log(row[0]);

                if (row[0] == selected_barcode) {

                    data_barcode.push(row);
                }
            });

            indexes.forEach(function(row) {
                if (row[0] == selected_barcode) {
                    //console.log(row);
                    category_labels.push(row);
                }

            });

            data_barcode.sort(function(a, b) {
                return a[2] - b[2];
            });

            data_barcode.forEach(function(row) {
                if (row[1].endsWith("Pass")){
                    //console.log(row);
                    categories.push(row);
                }


            })
            //console.log(data_barcode[-1]);
            //console.log(categories.map(x => x[2]));
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
                    categories: categories
                }
            });

            for (var i = 0; i < result_read_count_sum.length; i++) {

                    chart.addSeries({
                        name: result_read_count_sum[i]['name'],
                        data: result_read_count_sum[i]['data']
                    });
                //}
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
                    categories: categories
                }
            });

        for (var i = 0; i < result_read_length_sum.length; i++) {


                chart.addSeries({
                    name: result_read_length_sum[i]['name'],
                    data: result_read_length_sum[i]['data']
                });
        }
        //}
    });
}
