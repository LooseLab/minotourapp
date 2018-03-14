function requestChannelSummaryData(id) {
    /*
     * Request channel summary data
     */
    var url = "/api/v1/flowcells/" + id + "/channelsummary_readcount";

    $.get(url, function (data) {

        if (data.length > 0) {

            if (!self.chartReadsPerPore) {

                self.chartReadsPerPore = self.makeHeatmapChart(
                    "reads-per-pore",
                    "Reads per Channel".toUpperCase(),
                    ""
                );

            }

            self.updatePoreChart(self.chartReadsPerPore, data, 'read_count');

        }

    });

    var url = "/api/v1/flowcells/" + id + "/channelsummary_readkb";

    $.get(url, function (data) {

        if (data.length > 0) {

            if (!self.chartBasesPerPore) {

                self.chartBasesPerPore = self.makeHeatmapChart(
                    "bases-per-pore",
                    "bases (kb) per Channel".toUpperCase(),
                    ""
                );

            }

            self.updatePoreChart(self.chartBasesPerPore, data, 'read_length');

        }

    });

};
