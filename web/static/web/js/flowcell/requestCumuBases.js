function updateCumuBaseChart() {
};

function requestCumuBases(id) {

    var url = "/api/v1/flowcells/" + id + "/summarybarcodebyminute_bases";

    $.get(url, function (data) {

        if (!self.chart_cumulative_yield_overtime_new) {

            self.chart_cumulative_yield_overtime_new = self.makeChart2(
                "cumulative-yield-overtime-new",
                "cumulative bases".toUpperCase(),
                "cumulative bases".toUpperCase()
            );

        }

        var chart = self.chart_cumulative_yield_overtime_new;

        if (!self.chart_cumulative_number_reads_overtime_new) {

            self.chart_cumulative_number_reads_overtime_new = self.makeChart2(
                "cumulative-number-reads-overtime-new",
                "cumulative reads".toUpperCase(),
                "cumulative reads".toUpperCase()
            );

        }

        var chart2 = self.chart_cumulative_number_reads_overtime_new;

        var selectedBarcode = self.getSelectedBarcode();

        if (chart) {

            while (chart.series.length > 0) {
                chart.series[0].remove();
            }

        }

        if (chart2) {

            while (chart2.series.length > 0) {
                chart2.series[0].remove();
            }

        }

        for (var typeName of Object.keys(data[selectedBarcode]["bases"])) {
            for (var status of Object.keys(data[selectedBarcode]["bases"][typeName])) {
                chart.addSeries({
                    name: selectedBarcode + " - " + typeName + " - " + status,
                    data: data[selectedBarcode]["bases"][typeName][status]
                });
            }
        }

        for (var typeName of Object.keys(data[selectedBarcode]["reads"])) {
            for (var status of Object.keys(data[selectedBarcode]["reads"][typeName])) {
                chart2.addSeries({
                    name: selectedBarcode + " - " + typeName + " - " + status,
                    data: data[selectedBarcode]["reads"][typeName][status]
                });
            }
        }

        for (var i in self.rundata) {

            var starttime = new Date(Date.parse(self.rundata[i]['start_time']));
            var endtime = new Date(Date.parse(self.rundata[i]['last_read']));
            var name = self.rundata[i]['id']

            chart.xAxis[0].addPlotLine({
                value: starttime,
                color: 'black',
                dashStyle: 'dot',
                width: 2,
            })

            chart2.xAxis[0].addPlotLine({
                value: starttime,
                color: 'black',
                dashStyle: 'dot',
                width: 2,
            })

        }

    })

};
