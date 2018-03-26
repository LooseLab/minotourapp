function requestMaxLengthTimeCallback(data) {

    if (!this.max_read_lengths_overtime_new) {

        this.max_read_lengths_overtime_new = this.makeChart2(
            "max-read-lengths-overtime-new",
            "max read length over time".toUpperCase(),
            "max read length".toUpperCase()
        );

    }

    var chart = this.max_read_lengths_overtime_new;

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

function requestMaxLengthTime(id) {

    var url = "/api/v1/flowcells/" + id + "/summarybarcodebyminute_maxlength";

    requestMaxLengthTimeCallback = requestMaxLengthTimeCallback.bind(this);

    $.get(url, requestMaxLengthTimeCallback);

};

