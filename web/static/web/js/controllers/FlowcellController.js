class FlowcellController {
/*
This controller is the top level controller, initialised when the flowcell index page is loaded,
it contains references to the other controllers
 */
    constructor(flowcellId) {
        // delete soon
        this._interval = setInterval(requestData, 60000, flowcellId);
        this._flowcellTabController = new FlowcellTabController(flowcellId);
        this._lastread = 0; // TODO this is a quick fix for the problem on the live data tab. we should refactor this soon.
    }

    get flowcellTabController() {
        return this._flowcellTabController;
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
}
