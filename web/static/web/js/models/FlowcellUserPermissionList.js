class FlowcellUserPermissionList {
    /* Contains the list of flowcell user permissions */

    constructor() {
        
        this._flowcelluserpermissions = [];
    }

    adiciona(flowcelluserpermission) {
        
        this._flowcelluserpermissions.push(flowcelluserpermission);
    }

    get flowcelluserpermissions() {
        
        return [].concat(this._flowcelluserpermissions);
    }
}
