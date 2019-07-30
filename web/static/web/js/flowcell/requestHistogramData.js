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
                var chart_read_length = this.chartHistogramBasesSequencedByReadLength;

                if (!this.collectchartHistogramBasesSequencedByReadLength) {
                    this.collectchartHistogramBasesSequencedByReadLength = this.makeSplineChart(
                        "collect-histogram-bases-sequenced-by-read-length",
                        "Collect Histogram of Bases Sequenced by Read Length".toUpperCase(),
                        "Proportion of bases longer than x".toUpperCase()
                    );
                }

                var collect_chart_read_length = this.collectchartHistogramBasesSequencedByReadLength;

                if (!this.chartHistogramReadLength) {
                    this.chartHistogramReadLength = this.makeSplineChart(
                        "histogram-read-lengths",
                        "Histogram of Read Lengths".toUpperCase(),
                        "Number of Reads".toUpperCase()
                    );
                }
                var chart_read_count = this.chartHistogramReadLength;

                if (!this.collectchartHistogramReadLength) {
                    this.collectchartHistogramReadLength = this.makeSplineChart(
                        "collect-histogram-read-lengths",
                        "Collect Histogram of Read Lengths".toUpperCase(),
                        "Proportion of data longer than x".toUpperCase()
                    );
                }

                var collect_chart_read_count = this.collectchartHistogramReadLength;

                if (chart_read_count.series) {
                    while (chart_read_count.series.length > 0) {
                        chart_read_count.series[0].remove();
                    }
                }

                if (collect_chart_read_count.series) {
                    while (collect_chart_read_count.series.length > 0) {
                        collect_chart_read_count.series[0].remove();
                    }
                }

                if (chart_read_length.series) {
                    while (chart_read_length.series.length > 0) {
                        chart_read_length.series[0].remove();
                    }
                }

                if (collect_chart_read_length.series) {
                    while (collect_chart_read_length.series.length > 0) {
                        collect_chart_read_length.series[0].remove();
                    }
                }


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

                var options2 = {
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
                    },
                    yAxis: {
                        max: 100
                    }
                }

                chart_read_count.update(options);
                collect_chart_read_count.update(options2);
                chart_read_length.update(options);
                collect_chart_read_length.update(options2);

                Object.keys(data).forEach(function (barcode_name) {

                    Object.keys(data[barcode_name]).forEach(function(read_type_name){

                        Object.keys(data[barcode_name][read_type_name]).forEach(function (rejection_status) {

                            Object.keys(data[barcode_name][read_type_name][rejection_status]).forEach(function (status) {

                                chart_read_count.addSeries(
                                    data[barcode_name][read_type_name][rejection_status][status]["read_count"]
                                );

                                collect_chart_read_count.addSeries(
                                    data[barcode_name][read_type_name][rejection_status][status]["collect_read_count"]
                                );

                                chart_read_length.addSeries(
                                    data[barcode_name][read_type_name][rejection_status][status]["read_length"]
                                );

                                collect_chart_read_length.addSeries(
                                    data[barcode_name][read_type_name][rejection_status][status]["collect_read_length"]
                                );
                            });
                        });
                    });
                });
            }

        ).bind(this)
    )
    ;
}
