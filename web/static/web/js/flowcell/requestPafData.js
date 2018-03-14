
function requestPafData(id) {

    self.requestMappedChromosomes(id);

    var pafurl = '/api/v1/flowcells/' + id + '/pafsummary/';

    $.get(pafurl, function (data) {

        if (data.length > 0) {

            summarycoverage = {};

            for (var i = 0; i < data.length; i++) {

                if (summarycoverage[data[i].barcode_group_name] === undefined) {
                    summarycoverage[data[i].barcode_group_name] = {};
                }

                if (summarycoverage[data[i].barcode_group_name][data[i].read_type_name] === undefined) {
                    summarycoverage[data[i].barcode_group_name][data[i].read_type_name] = {};
                }

                if (summarycoverage[data[i].barcode_group_name][data[i].read_type_name][data[i].chrom_name] === undefined) {
                    summarycoverage[data[i].barcode_group_name][data[i].read_type_name][data[i].chrom_name] = {};
                }

                summarycoverage[data[i].barcode_group_name][data[i].read_type_name][data[i].chrom_name]["coverage"] = {
                    "name": "coverage",
                    "data": [data[i].chrom_cover],
                    "animation": false
                };

                summarycoverage[data[i].barcode_group_name][data[i].read_type_name][data[i].chrom_name]["ave_read_len"] = {
                    "name": "Average Read Length",
                    "data": [data[i].avg_read_len],
                    "animation": false
                };

            }

            self.summarycoverage = summarycoverage;
            self.updateCoverageBasedCharts(self.chart_per_chrom_cov, "coverage");
            self.updateCoverageBasedCharts(self.chart_per_chrom_avg, "ave_read_len");

        }

    })

}
