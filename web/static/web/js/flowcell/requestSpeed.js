function requestSpeed(id) {

    var url = "/api/v1/flowcells/" + id + "/summarybarcodebyminute_speed";

    var self = this;

    var selectedBarcode = this.getSelectedBarcode();

    $.get(url, function (data) {



        if (!self.chartSequencingRate_new) {
            self.chartSequencingRate_new = self.makeChart2(
                "sequencing-rate-new",
                "sequencing rate".toUpperCase(),
                "bases/second".toUpperCase()
            );
        }

        var chart = self.chartSequencingRate_new;

        if (!self.chartSequencingSpeed_new) {
            self.chartSequencingSpeed_new = self.makeChart2(
                "sequencing-speed-new",
                "sequencing speed".toUpperCase(),
                "bases/channel/second".toUpperCase()
            );
        }

        var chart2 = self.chartSequencingSpeed_new;

        if (chart.series) {
            while (chart.series.length > 0) {
                chart.series[0].remove();
            }
        }

        if (chart2.series) {
            while (chart2.series.length > 0) {
                chart2.series[0].remove();
            }
        }

        for (var typeName of Object.keys(data[selectedBarcode]["speed"])) {
            for (var status of Object.keys(data[selectedBarcode]["speed"][typeName])) {
                chart2.addSeries({
                    name: selectedBarcode + " - " + typeName + " - " + status,
                    data: data[selectedBarcode]["speed"][typeName][status]
                });
            }
        }

        for (var typeName of Object.keys(data[selectedBarcode]["rate"])) {
            for (var status of Object.keys(data[selectedBarcode]["rate"][typeName])) {
                chart.addSeries({
                    name: selectedBarcode + " - " + typeName + " - " + status,
                    data: data[selectedBarcode]["rate"][typeName][status]
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
