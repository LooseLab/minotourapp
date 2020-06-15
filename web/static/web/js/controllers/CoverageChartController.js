class CoverageChartController {
    // controller for coverage chart
    constructor(div_name) {
        // constructs new Coverage chart class, and link to chromosome select field
        this._chromosomeSelect = document.querySelector('#chromosome-id-select');
        this._coverageChart = new CoverageChart(div_name);
    }

    get coverage_chart() {
        return this._coverageChart;
    }

    reload_master_chart() {
        // reload the data on click to reset or chromosome change
        let url = this.create_url();
        this.loadChartData(url);
    }

    create_url() {

        // create the url to query for data
        var selected_option = this._chromosomeSelect.value;

        if (parseInt(selected_option) < 0) {

            return "/api/v1/pafcoverage/0/0/0/1/";
        }

        // get the selected option to view coverage
        var value_combination = selected_option.split('_');
        console.log(value_combination);
        // get the task id for the jobmaster
        var task_id = value_combination[0];
        // get the barcode name
        var barcode_name = value_combination[1];
        // template
        var read_type_id = value_combination[2];
        // get the chromosome name
        var chromosome_id = value_combination[3];
        // create the url
        var url = "/api/v1/pafcoverage/" + task_id + "/" + barcode_name + "/" + read_type_id + "/" + chromosome_id + "/";

        return url;
    }

    loadChartData(url) {

        let self = this;

        $.getJSON(url, (function (data) {
            var parsedData = JSON.parse(data);
            self._coverageChart.masterChart.series[0].setData(parsedData);
            self._coverageChart.detailChart.series[0].setData(parsedData);
            self._coverageChart.masterChart.xAxis[0].removePlotBand('mask-before');
            return data;
        }));
    }
}