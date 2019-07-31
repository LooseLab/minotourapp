class CoverageChartController {
    // controller for coverage chart
    constructor(div_name, advanced) {
        // constructs new Coverage chart class, and link to chromosome select field
        this._chromosome_select = document.querySelector('#' + advanced + 'chromosome-id-select');
        if(advanced === ""){
            this._coverage_chart = new CoverageChart(div_name);
        }
    }

    get coverage_chart() {
        return this._coverage_chart;
    }

    reload_master_chart() {
        // reload the data on click to reset or chromosome change
        let url = this.create_url();
        this.load_chart_data(url);
    }

    create_url() {

        // create the url to query for data
        var selected_option = this._chromosome_select.value;

        if (parseInt(selected_option) < 0) {

            return "/api/v1/pafcoverage/0/0/0/1/";
        }

        // get the selected option to view coverage
        var value_combination = selected_option.split('_');
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

    load_chart_data(url) {

        let self = this;

        $.getJSON(url, (function (data) {
            var parsedData = JSON.parse(data);
            self._coverage_chart.master_chart.series[0].setData(parsedData);
            self._coverage_chart.detail_chart.series[0].setData(parsedData);
            self._coverage_chart.master_chart.xAxis[0].removePlotBand('mask-before');
            return data;
        }));
    }
}