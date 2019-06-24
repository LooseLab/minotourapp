class FlowcellService {

    constructor() {
        this._http = new HttpService();
    }

    getFlowcellTabs(flowcell_id) {

        return this._http
            .get(`/api/v1/flowcells/${flowcell_id}/tabs/`)
            .then(tabs => {

                return tabs;
            })
            .catch(erro => {
                console.log(erro);
                throw new Error('Error getting the flowcell tabs.');
            });
    }
}