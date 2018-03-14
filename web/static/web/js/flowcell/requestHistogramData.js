function requestHistogramData(id) {
    /*
     * Request histogram data
     */

    var url = "/api/v1/flowcells/" + id + "/histogramsummary";

    $.get(url, function (data) {
        //console.log(data);
        if (data.length > 0) {

            var ordered_data = data.sort(function (a, b) {
                return a.bin_width - b.bin_width;
            });

            var summary = {};
            //console.log(ordered_data);
            for (var i = 0; i < ordered_data.length; i++) {

                var item = ordered_data[i];
                //console.log(item);
                if (summary[item.barcode_name] === undefined) {
                    summary[item.barcode_name] = {};
                }

                if (summary[item.barcode_name][item.read_type_name] === undefined) {
                    summary[item.barcode_name][item.read_type_name] = {};
                }

                if (summary[item.barcode_name][item.read_type_name][parseInt(item.bin_width)] === undefined) {
                    summary[item.barcode_name][item.read_type_name][parseInt(item.bin_width)] = {};
                    summary[item.barcode_name][item.read_type_name][parseInt(item.bin_width)]['read_count'] = 0;
                    summary[item.barcode_name][item.read_type_name][parseInt(item.bin_width)]['read_length'] = 0;
                }
                summary[item.barcode_name][item.read_type_name][parseInt(item.bin_width)]['read_count'] += item.read_count;
                summary[item.barcode_name][item.read_type_name][parseInt(item.bin_width)]['read_length'] += item.read_length;
            }
            ;
            //console.log("summary");
            //console.log(summary);
            var summaries = {};

            for (var barcode in summary) {

                for (var read_type in summary[barcode]) {
                    var lastbinwidth = 0;

                    //Here we make the empty arrays
                    if (summaries[barcode] === undefined) {
                        summaries[barcode] = {};
                    }
                    if (summaries[barcode][read_type] === undefined) {
                        summaries[barcode][read_type] = {
                            'bin_width': [],
                            'read_count': [],
                            'read_length': [],
                        };
                    }

                    var maxbinwidth = 0;

                    for (var tmp_bin_width of Object.keys(summary[barcode][read_type])) {
                        if (parseInt(tmp_bin_width) > maxbinwidth) {
                            maxbinwidth = parseInt(tmp_bin_width);
                        }
                    }

                    //var maxbinwidth = Math.max(Object.keys(summary[barcode][read_type]))

                    //var maxbinwidth = summary[barcode][read_type]
                    for (var i = 0; i <= maxbinwidth; i++) {
                        summaries[barcode][read_type]['bin_width'].push(parseInt(i * 900 + 900));
                        summaries[barcode][read_type]['read_count'].push(0);
                        summaries[barcode][read_type]['read_length'].push(0);
                    }
                    for (var bin_width in summary[barcode][read_type]) {
                        summaries[barcode][read_type]['read_count'][bin_width] = summary[barcode][read_type][bin_width]['read_count'];
                        summaries[barcode][read_type]['read_length'][bin_width] = summary[barcode][read_type][bin_width]['read_length'];
                        summaries[barcode][read_type]['bin_width'][bin_width] = parseInt(bin_width * 900 + 900);

                    }
                }
            }


            self.histogramSummary = summaries;

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

            for (var barcode_name of Object.keys(summaries)) {

                if (barcode_name === self.selectedBarcode) {

                    for (var typeName of Object.keys(summaries[barcode_name])) {

                        chart.update({
                            chart: {
                                type: 'column'
                            },
                            xAxis: {
                                type: 'category',
                                categories: summaries[barcode_name][typeName]['bin_width']
                            }
                        });

                        chart.addSeries({
                            name: barcode_name + " - " + typeName,
                            data: summaries[barcode_name][typeName]["read_length"]
                        });

                    }
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

            for (var barcode_name of Object.keys(summaries)) {

                if (barcode_name === self.selectedBarcode) {

                    for (var typeName of Object.keys(summaries[barcode_name])) {

                        chart.update({
                            chart: {
                                type: 'column'
                            },
                            xAxis: {
                                type: 'category',
                                categories: summaries[barcode_name][typeName]['bin_width']
                            }
                        });

                        chart.addSeries({
                            name: barcode_name + " - " + typeName,
                            data: summaries[barcode_name][typeName]["read_count"]
                        });

                    }
                }
            }


        }

    });

}
