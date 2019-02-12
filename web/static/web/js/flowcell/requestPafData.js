function updateCoverageBasedCharts(chart, summarycoverage, field) {

    var series = [];
    var categories = [];
    for (var barcode of Object.keys(summarycoverage)) {

        data = [];
        for (var readtype of Object.keys(summarycoverage[barcode])) {

            for (var chromosome of Object.keys(summarycoverage[barcode][readtype])) {

                categories.push(chromosome);
                data.push(summarycoverage[barcode][readtype][chromosome][field]["data"]);
            }
        }
        serie = {
            "name": barcode + " " + readtype,
            "data": data
        };
        series.push(serie);
    }

    chart.xAxis[0].setCategories(categories);
    var chartSeriesLength = (chart.series ? chart.series.length : 0);
    for (var i = 0; i < series.length; i++) {
        if (i <= (chartSeriesLength - 1)) {
            chart.series[i].setData(series[i].data);
            chart.series[i].update({
                name: series[i].name
            });
        } else {
            chart.addSeries(series[i]);
        }
    }
}

function requestMappedChromosomes(flowcell_id) {
    /*
     * Request the chromosomes that have reads mapped to using minimap2
     * and update the select box on tab Mapping
     *
     */

    var url = '/api/v1/flowcells/' + flowcell_id + '/references';

    $.getJSON(url, function (data_unsorted) {

        var data = data_unsorted.sort(function (a, b) {
            return a['barcode_name'] - ['b.barcode_name'];
        });

        var select = document.getElementById('chromosome-id-select');
        var selected_index = select.selectedIndex;
        var selected_option = select[selected_index];

        while (select.length > 0) {
            select.remove(0);
        }

        var option = document.createElement('option');
        option.text = '--- Select ---';
        option.value = '-1';
        select.add(option);

        for (var i = 0; i < data.length; i++) {
            var option = document.createElement('option');

            var value_combination = data[i]['task_id'] + '_' + data[i]['barcode_name'] + '_' + data[i]['read_type_id'] + '_' + data[i]['chromosome_id'];

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

function requestPafData(id) {

    requestMappedChromosomes(id);

    var pafurl = '/api/v1/flowcells/' + id + '/pafsummary/';

    if (!this.chart_per_chrom_cov) {

        this.chart_per_chrom_cov = this.makeChart(
            "per-chrom-cov",
            "Chromosome Coverage".toUpperCase(),
            "Chromosome Coverage".toUpperCase()
        );

    }

    if (!this.chart_per_chrom_avg) {

        this.chart_per_chrom_avg = this.makeChart(
            "per-chrom-avg",
            "Read Length By Chromosome".toUpperCase(),
            "Read Length By Chromosome".toUpperCase()
        );

    }

    $.get(pafurl, (function (data_obj) {

        var data = data_obj.data;

        if (data.length > 0) {

            summarycoverage = {};

            for (var i = 0; i < data.length; i++) {

                if (summarycoverage[data[i].barcode_name] === undefined) {
                    summarycoverage[data[i].barcode_name] = {};
                }

                if (summarycoverage[data[i].barcode_name][data[i].read_type_name] === undefined) {
                    summarycoverage[data[i].barcode_name][data[i].read_type_name] = {};
                }

                if (summarycoverage[data[i].barcode_name][data[i].read_type_name][data[i].reference_line_name] === undefined) {
                    summarycoverage[data[i].barcode_name][data[i].read_type_name][data[i].reference_line_name] = {};
                }

                summarycoverage[data[i].barcode_name][data[i].read_type_name][data[i].reference_line_name]["coverage"] = {
                    "name": "coverage",
                    "data": [data[i].chrom_cover],
                    "animation": false
                };

                summarycoverage[data[i].barcode_name][data[i].read_type_name][data[i].reference_line_name]["ave_read_len"] = {
                    "name": "Average Read Length",
                    "data": [data[i].avg_read_len],
                    "animation": false
                };

            }

            updateCoverageBasedCharts(this.chart_per_chrom_cov, summarycoverage, "coverage");
            updateCoverageBasedCharts(this.chart_per_chrom_avg, summarycoverage, "ave_read_len");

        }

    }).bind(this));

    /*
     * Request flowcell summary
     * This response is a HTML code that is appended to the page (no json)
     */
    var url = "/flowcells/" + id + "/pafsummaryhtml";

    var result_mapping_summary = document.querySelector('#mapping-summary');

    $.ajax({
        url: url,
        dataType: "html",
        success: function(data) {
            result_mapping_summary.innerHTML = data;
        },
        error: function(e)
        {
            console.log('Error: ' + e);
        }
    });

}
