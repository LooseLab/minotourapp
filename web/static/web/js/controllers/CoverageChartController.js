class CoverageChartController {

    constructor(div_name) {

        this._chromossome_select = document.querySelector('#chromosome-id-select');
        this._coverage_chart = new CoverageChart(div_name);
    }

    get coverage_chart() {

        return this._coverage_chart;
    }

    reload_master_chart() {

        let url = this.create_url();
        this.load_chart_data(url);
    }

    create_url() {

        // var selected_index = this.select_container.selectedIndex;
        var selected_option = this._chromossome_select.value;

        if (parseInt(selected_option) < 0) {

            return "/api/v1/pafcoverage/0/0/0/1/";
        }

        // var value_combination = data[i]['task_id'] + '_' + data[i]['barcode_name'] + '_' + data[i]['read_type_id'] + '_' + '_' + data[i]['chromosome_id'];

        var value_combination = selected_option.split('_');
        var task_id = value_combination[0];
        // var run_id = value_combination[0];
        var barcode_name = value_combination[1];
        // var barcode_id = value_combination[1];
        var read_type_id = value_combination[2];
        // var reference_id = value_combination[3];
        var chromosome_id = value_combination[3];
        // var chromosome_id = value_combination[4];

        // var url = "/api/v1/flowcells/pafcoverage/" + task_id + "/" + barcode_name + "/" + read_type_id + "/" + chromosome_id + "/";
        var url = "/api/v1/pafcoverage/" + task_id + "/" + barcode_name + "/" + read_type_id + "/" + chromosome_id + "/";

        return url;
    }

    load_chart_data(url) {

        let self = this;

        $.getJSON(url, (function (datax) {

            var data = JSON.parse(datax);

            self._coverage_chart.master_chart.series[0].setData(data);
            self._coverage_chart.detail_chart.series[0].setData(data);
            self._coverage_chart.master_chart.xAxis[0].removePlotBand('mask-before');

            return data;
        }));
    }
}