function updateCumuBaseChart() {
};

function requestCumuBasesCallBack(data) {

    if (!this.chart_cumulative_yield_overtime_new) {

        this.chart_cumulative_yield_overtime_new = this.makeChart2(
            "cumulative-yield-overtime-new",
            "cumulative bases".toUpperCase(),
            "cumulative bases".toUpperCase()
        );

    }

    var chart = this.chart_cumulative_yield_overtime_new;

    if (!this.chart_cumulative_number_reads_overtime_new) {

        this.chart_cumulative_number_reads_overtime_new = this.makeChart2(
            "cumulative-number-reads-overtime-new",
            "cumulative reads".toUpperCase(),
            "cumulative reads".toUpperCase()
        );

    }

    var chart2 = this.chart_cumulative_number_reads_overtime_new;

    var selected_barcode = get_selected_barcode();

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

    for (var typeName of Object.keys(data[selected_barcode]["bases"])) {

        for (var status of Object.keys(data[selected_barcode]["bases"][typeName])) {

            chart.addSeries({
                name: selected_barcode + " - " + typeName + " - " + status,
                data: data[selected_barcode]["bases"][typeName][status]
            });
        }
    }

    for (var typeName of Object.keys(data[selected_barcode]["reads"])) {
        for (var status of Object.keys(data[selected_barcode]["reads"][typeName])) {
            chart2.addSeries({
                name: selected_barcode + " - " + typeName + " - " + status,
                data: data[selected_barcode]["reads"][typeName][status]
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

        chart2.xAxis[0].addPlotLine({
            value: starttime,
            color: 'black',
            dashStyle: 'dot',
            width: 2,
        })

    }

}


function requestCumuBases(id) {

    var url = "/api/v1/flowcells/" + id + "/summarybarcodebyminute_bases";

    requestCumuBasesCallBack = requestCumuBasesCallBack.bind(this);

    $.get(url, requestCumuBasesCallBack);

};
