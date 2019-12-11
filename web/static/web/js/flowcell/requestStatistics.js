function requestStatisticsCallback(data) {
    // if we don't have any data
    if (!data && data.length <= 0) {
        return;
    }
    // if we don't already have this chart, initialise it
    if (!this.average_quality_overtime_new) {

        this.average_quality_overtime_new = this.makeSplineChart(
            "average-quality-overtime-new",
            "average quality over time".toUpperCase(),
            "average read quality score".toUpperCase()
        );

    }
// if we don't already have this chart, initialise it
    if (!this.average_read_lengths_overtime_new) {

        this.average_read_lengths_overtime_new = this.makeSplineChart(
            "average-read-lengths-overtime-new",
            "average read length over time".toUpperCase(),
            "average read length".toUpperCase()
        );

    }
// if we don't already have this chart, initialise it
    if (!this.chart_cumulative_yield_overtime_new) {

        this.chart_cumulative_yield_overtime_new = this.makeSplineChart(
            "cumulative-yield-overtime-new",
            "cumulative bases".toUpperCase(),
            "cumulative bases".toUpperCase()
        );

    }
// if we don't already have this chart, initialise it
    if (!this.chart_cumulative_number_reads_overtime_new) {

        this.chart_cumulative_number_reads_overtime_new = this.makeSplineChart(
            "cumulative-number-reads-overtime-new",
            "cumulative reads".toUpperCase(),
            "cumulative reads".toUpperCase()
        );

    }
// if we don't already have this chart, initialise it
    if (!this.max_read_lengths_overtime_new) {

        this.max_read_lengths_overtime_new = this.makeSplineChart(
            "max-read-lengths-overtime-new",
            "max read length over time".toUpperCase(),
            "max read length".toUpperCase()
        );

    }
// if we don't already have this chart, initialise it
    if (!this.chartSequencingRate_new) {
        this.chartSequencingRate_new = this.makeSplineChart(
            "sequencing-rate-new",
            "sequencing rate".toUpperCase(),
            "bases/second".toUpperCase()
        );
    }


    var average_quality_chart = this.average_quality_overtime_new;

    var average_read_length_chart = this.average_read_lengths_overtime_new;

    var cumulative_yield_chart = this.chart_cumulative_yield_overtime_new;

    var cumulative_reads_chart = this.chart_cumulative_number_reads_overtime_new;

    var max_read_length_chart = this.max_read_lengths_overtime_new;

    var sequencing_rate_chart = this.chartSequencingRate_new;

    //var chart7 = this.chartSequencingSpeed_new;

    // var charts = [chart, chart2, chart3, chart4, chart5, sequencing_rate_chart, chart7];
    var charts = [average_quality_chart, average_read_length_chart, cumulative_yield_chart,
        cumulative_reads_chart, max_read_length_chart, sequencing_rate_chart];

    var selected_barcode = get_selected_barcode();

    // for each chart, remove any existing data
    charts.forEach(function(value) {

        if (value) {

            while (value.series.length > 0) {
                value.series[0].remove();
            }
        }
    });

    var rundata = data['runs'];

    var data_keys = data['data_keys'];

    var i = 0;
    for (i = 0; i < data_keys.length; i++) {
        var data_average_quality_over_time = data['data'][data_keys[i]].map(x => [x[0], x[1]]);
        var data_average_read_length_over_time = data['data'][data_keys[i]].map(x => [x[0], x[2]]);
        var data_cumulative_bases = data['data'][data_keys[i]].map(x => [x[0], x[3]]);
        var data_cumulative_reads = data['data'][data_keys[i]].map(x => [x[0], x[4]]);
        var data_max_length = data['data'][data_keys[i]].map(x => [x[0], x[5]]);
        var data_sequencing_rate = data['data'][data_keys[i]].map(x => [x[0], x[6]]);
        //var data_sequencing_speed = data['data'][data_keys[i]].map(x => [x[0], x[7]]);

        average_quality_chart.addSeries({
            name: data_keys[i],
            data: data_average_quality_over_time
        });

        average_read_length_chart.addSeries({
            name: data_keys[i],
            data: data_average_read_length_over_time
        });

        cumulative_yield_chart.addSeries({
            name: data_keys[i],
            data: data_cumulative_bases
        });

        cumulative_reads_chart.addSeries({
            name: data_keys[i],
            data: data_cumulative_reads
        });

        max_read_length_chart.addSeries({
            name: data_keys[i],
            data: data_max_length
        });

        sequencing_rate_chart.addSeries({
            name: data_keys[i],
            data: data_sequencing_rate
        });

        //chart7.addSeries({
        //    name: data_keys[i],
        //    data: data_sequencing_speed
        //});

    }

    for (var i in rundata) {

        var starttime = new Date(Date.parse(rundata[i]['start_time']));
        var endtime = new Date(Date.parse(rundata[i]['last_read']));
        var name = rundata[i]['id'];

        average_quality_chart.xAxis[0].addPlotLine({
            value: starttime,
            color: 'black',
            dashStyle: 'dot',
            width: 2,
        });

        average_read_length_chart.xAxis[0].addPlotLine({
            value: starttime,
            color: 'black',
            dashStyle: 'dot',
            width: 2,
        });

        cumulative_yield_chart.xAxis[0].addPlotLine({
            value: starttime,
            color: 'black',
            dashStyle: 'dot',
            width: 2,
        });

        cumulative_reads_chart.xAxis[0].addPlotLine({
            value: starttime,
            color: 'black',
            dashStyle: 'dot',
            width: 2,
        });

        max_read_length_chart.xAxis[0].addPlotLine({
            value: starttime,
            color: 'black',
            dashStyle: 'dot',
            width: 2,
        });

        sequencing_rate_chart.xAxis[0].addPlotLine({
            value: starttime,
            color: 'black',
            dashStyle: 'dot',
            width: 2,
        });

        //chart7.xAxis[0].addPlotLine({
        //    value: starttime,
        //    color: 'black',
        //    dashStyle: 'dot',
        //    width: 2,
        //});
    }
}

function requestStatistics(id) {

    var selected_barcode = get_selected_barcode();

    var url = "/api/v1/flowcells/" + id + "/statistics/?barcode_name=" + selected_barcode;;

    requestStatisticsCallback = requestStatisticsCallback.bind(this);

    $.get(url, requestStatisticsCallback);

};
