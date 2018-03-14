function requestMaxLengthTime(id) {

    var url = "/api/v1/flowcells/" + id + "/summarybarcodebyminute_maxlength";

    $.get(url, function (data) {

        if (!self.max_read_lengths_overtime_new) {

            self.max_read_lengths_overtime_new = self.makeChart2(
                "max-read-lengths-overtime-new",
                "max read length over time".toUpperCase(),
                "max read length".toUpperCase()
            );

        }

        var chart = self.max_read_lengths_overtime_new;

        var selectedBarcode = self.getSelectedBarcode();

        if (chart) {

            while (chart.series.length > 0) {
                chart.series[0].remove();
            }

        }

        for (var typeName of Object.keys(data[selectedBarcode])) {
            for (var status of Object.keys(data[selectedBarcode][typeName])) {
                chart.addSeries({
                    name: selectedBarcode + " - " + typeName + " - " + status,
                    data: data[selectedBarcode][typeName][status]
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
        }

    })

};

