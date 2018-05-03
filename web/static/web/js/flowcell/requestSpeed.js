function requestSpeedCallback(data) {

    var selected_barcode = get_selected_barcode();

    if (!this.chartSequencingRate_new) {
        this.chartSequencingRate_new = this.makeChart2(
            "sequencing-rate-new",
            "sequencing rate".toUpperCase(),
            "bases/second".toUpperCase()
        );
    }

    var chart = this.chartSequencingRate_new;

    if (!this.chartSequencingSpeed_new) {
        this.chartSequencingSpeed_new = this.makeChart2(
            "sequencing-speed-new",
            "sequencing speed".toUpperCase(),
            "bases/channel/second".toUpperCase()
        );
    }

    var chart2 = this.chartSequencingSpeed_new;

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

    for (var typeName of Object.keys(data['speed'][selected_barcode])) {
        for (var status of Object.keys(data['speed'][selected_barcode][typeName])) {
            chart2.addSeries({
                name: selected_barcode + " - " + typeName + " - " + status,
                data: data['speed'][selected_barcode][typeName][status]
            });
        }
    }

    for (var typeName of Object.keys(data["rate"][selected_barcode])) {
        for (var status of Object.keys(data["rate"][selected_barcode][typeName])) {
            chart.addSeries({
                name: selected_barcode + " - " + typeName + " - " + status,
                data: data["rate"][selected_barcode][typeName][status]
            });
        }
    }

    // for (var i in self.rundata) {
    //     var starttime = new Date(Date.parse(self.rundata[i]['start_time']));
    //     var endtime = new Date(Date.parse(self.rundata[i]['last_read']));
    //     var name = self.rundata[i]['id']
    //     chart.xAxis[0].addPlotLine({
    //         value: starttime,
    //         color: 'black',
    //         dashStyle: 'dot',
    //         width: 2,
    //     })
    //     chart2.xAxis[0].addPlotLine({
    //         value: starttime,
    //         color: 'black',
    //         dashStyle: 'dot',
    //         width: 2,
    //     })
    // }


};

function requestSpeed(id) {

    var url = "/api/v1/flowcells/" + id + "/summarybarcodebyminute_speed/";

    requestSpeedCallback = requestSpeedCallback.bind(this);

    $.get(url, requestSpeedCallback);

};

