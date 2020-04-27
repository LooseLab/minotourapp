class FlowcellController {
/*
This controller is the top level controller, initialised when the flowcell index page is loaded,
it contains references to the other controllers
 */
    constructor(flowcell_id) {

        console.log('Initialising FlowcellController.');
        this._interval = setInterval(requestData, 60000, flowcell_id);
        this._flowcell = new Flowcell(flowcell_id);
        this._articController = new ArticController(flowcell_id);
        this._flowcell_tab_controller = new FlowcellTabController(flowcell_id);
        this._coverage_chart_controller = new CoverageChartController('coverage_div', "");
        this._lastread = 0; // TODO this is a quick fix for the problem on the live data tab. we should refactor this soon.
    }

    get flowcell_tab_controller() {

        return this._flowcell_tab_controller;
    }

    get coverage_chart_controller() {

        return this._coverage_chart_controller;
    }

    get lastread() {

        return this._lastread;
    }

    set lastread(val) {

        this._lastread = val;
    }

    get articController() {
        return this._articController;
    }

}
