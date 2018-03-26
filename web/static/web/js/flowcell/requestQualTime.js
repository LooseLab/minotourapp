function requestQualTimeCallback(data) {

    if (!this.average_quality_overtime_new) {

        this.average_quality_overtime_new = this.makeChart2(
            "average-quality-overtime-new",
            "average quality over time".toUpperCase(),
            "average read quality score".toUpperCase()
        );

    }

    var chart = this.average_quality_overtime_new;

    var selectedBarcode = this.getSelectedBarcode();

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

    for (var i in this.rundata) {
        var starttime = new Date(Date.parse(this.rundata[i]['start_time']));
        var endtime = new Date(Date.parse(this.rundata[i]['last_read']));
        var name = this.rundata[i]['id']
        chart.xAxis[0].addPlotLine({
            value: starttime,
            color: 'black',
            dashStyle: 'dot',
            width: 2,
        })
    }

}

function requestQualTime(id) {

    var url = "/api/v1/flowcells/" + id + "/summarybarcodebyminute_quality";

    requestQualTimeCallback = requestQualTimeCallback.bind(this);

    $.get(url, requestQualTimeCallback);

};
