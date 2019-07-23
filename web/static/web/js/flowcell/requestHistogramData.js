function requestHistogramData(id) {
    /*
     * Request histogram data - Request data for the two histograms on the basecalled data tab
     */

    var selected_barcode = get_selected_barcode();

    var url = "/api/v1/flowcells/" + id + "/histogramsummary/?barcode_name=" + selected_barcode;
    // ajax request for the histogram data
    $.get(url, (function (dataObj) {

        var data = dataObj['data'];

        var categories = dataObj['categories'];
        // if we don't already have a chart, set one up calling the makeSplineChart function
        if (!this.chartHistogramBasesSequencedByReadLength) {

            this.chartHistogramBasesSequencedByReadLength = this.makeSplineChart(
                "histogram-bases-sequenced-by-read-length",
                "Histogram of Bases Sequenced by Read Length".toUpperCase(),
                "Number of bases".toUpperCase()
            );

        }
        // set the yield Histogram to a local scope rather than editing the instance floating around globally
        var chart_read_length = this.chartHistogramBasesSequencedByReadLength;
        // if we don't have a total number of reads histogram
        if (!this.chartHistogramReadLength) {
            // make one
            this.chartHistogramReadLength = this.makeSplineChart(
                "histogram-read-lengths",
                "Histogram of Read Lengths".toUpperCase(),
                "Number of reads".toUpperCase()
            );

        }
        // assign local scope to the total reads histogram
        var chart_read_count = this.chartHistogramReadLength;
        // if there is already data, remove it all on both charts
        if (chart_read_count.series) {
            while (chart_read_count.series.length > 0) {
                chart_read_count.series[0].remove();
            }
        }

        if (chart_read_length.series) {
            while (chart_read_length.series.length > 0) {
                chart_read_length.series[0].remove();
            }
        }
        // a object containg all the plot options
        var options = {
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
        };
        // update the charts with new options
        chart_read_count.update(options);
        chart_read_length.update(options);

        // loop through the barcodes in the results dictionary
        Object.keys(data).forEach(function (barcode_name) {
            // next level dictionary, loop through the keys (template, 1d^2 etc.)
            Object.keys(data[barcode_name]).forEach(function (read_type_name) {
                // loop through the next level of keys (sequence, unblock)
                Object.keys(data[barcode_name][read_type_name]).forEach(function (rejection_status) {
                    // loop through pass fail dictionary keys
                    Object.keys(data[barcode_name][read_type_name][rejection_status]).forEach(function (status) {
                        // add series for all the tree paths through the dictionary to the two charts
                        chart_read_count.addSeries(
                            data[barcode_name][read_type_name][rejection_status][status]["read_count"]
                        );

                        chart_read_length.addSeries(
                            data[barcode_name][read_type_name][rejection_status][status]["read_length"]
                        );
                    });
                });
            });
        });
    }).bind(this));
}
