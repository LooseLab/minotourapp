function requestQualTimeCallback(data) {

    if (!data && data.length <= 0) {
        return;
    }

    if (!this.average_quality_overtime_new) {

        this.average_quality_overtime_new = this.makeChart2(
            "average-quality-overtime-new",
            "average quality over time".toUpperCase(),
            "average read quality score".toUpperCase()
        );

    }

    if (!this.average_read_lengths_overtime_new) {

        this.average_read_lengths_overtime_new = this.makeChart2(
            "average-read-lengths-overtime-new",
            "average read length over time".toUpperCase(),
            "average read length".toUpperCase()
        );

    }

    if (!this.chart_cumulative_yield_overtime_new) {

        this.chart_cumulative_yield_overtime_new = this.makeChart2(
            "cumulative-yield-overtime-new",
            "cumulative bases".toUpperCase(),
            "cumulative bases".toUpperCase()
        );

    }

    if (!this.chart_cumulative_number_reads_overtime_new) {

        this.chart_cumulative_number_reads_overtime_new = this.makeChart2(
            "cumulative-number-reads-overtime-new",
            "cumulative reads".toUpperCase(),
            "cumulative reads".toUpperCase()
        );

    }

    if (!this.max_read_lengths_overtime_new) {

        this.max_read_lengths_overtime_new = this.makeChart2(
            "max-read-lengths-overtime-new",
            "max read length over time".toUpperCase(),
            "max read length".toUpperCase()
        );

    }

    var chart = this.average_quality_overtime_new;

    var chart2 = this.average_read_lengths_overtime_new;

    var chart3 = this.chart_cumulative_yield_overtime_new;

    var chart4 = this.chart_cumulative_number_reads_overtime_new;

    var chart5 = this.max_read_lengths_overtime_new;

    var charts = [chart, chart2, chart3, chart4, chart5];

    var selected_barcode = get_selected_barcode();

    charts.forEach(function(value) {

        if (value) {

            while (value.series.length > 0) {
                value.series[0].remove();
            }
        }
    });

    var rundata = data['runs'];
    var data_keys = data['data_keys'];
    // var data = data['data'];

    var i = 0;
    for (i = 0; i < data_keys.length; i++) {

        var data_average_quality_over_time = data['data'][data_keys[i]].map(x => [x[0], x[1]]);
        var data_average_read_length_over_time = data['data'][data_keys[i]].map(x => [x[0], x[2]]);
        var data_cumulative_bases = data['data'][data_keys[i]].map(x => [x[0], x[3]]);
        var data_cumulative_reads = data['data'][data_keys[i]].map(x => [x[0], x[4]]);
        var data_max_length = data['data'][data_keys[i]].map(x => [x[0], x[5]]);

        chart.addSeries({
            name: data_keys[i],
            data: data_average_quality_over_time
        });

        chart2.addSeries({
            name: data_keys[i],
            data: data_average_read_length_over_time
        });

        chart3.addSeries({
            name: data_keys[i],
            data: data_cumulative_bases
        });

        chart4.addSeries({
            name: data_keys[i],
            data: data_cumulative_reads
        });

        chart5.addSeries({
            name: data_keys[i],
            data: data_cumulative_reads
        });

    }

    for (var i in rundata) {
        var starttime = new Date(Date.parse(rundata[i]['start_time']));
        var endtime = new Date(Date.parse(rundata[i]['last_read']));
        var name = rundata[i]['id']
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
        chart3.xAxis[0].addPlotLine({
            value: starttime,
            color: 'black',
            dashStyle: 'dot',
            width: 2,
        })
        chart4.xAxis[0].addPlotLine({
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
