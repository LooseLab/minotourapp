class FlowcellController {

    constructor(flowcell_id) {
        this._interval = setInterval(requestData, 60000, flowcell_id);
        this._flowcell = new Flowcell(flowcell_id);

        this._flowcell_tab_controller = new FlowcellTabController(flowcell_id);
        this._coverage_chart_controller = new CoverageChartController('coverage_div');
    }

    get flowcell_tab_controller() {

        return this._flowcell_tab_controller;
    }

    get coverage_chart_controller() {

        return this._coverage_chart_controller;
    }
}