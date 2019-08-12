class FlowcellUserPermission {
    
    constructor(flowcell_id, user_id, permission) {
        this._flowcell_id = flowcell_id;
        this._user_id = user_id;
        this._permission = permission;
    }

    get flowcell_id() {

        return this._flowcell_id;
    }
    
    get user_id() {

        return this._user_id;
    }
    
    get permission() {

        return this.permission;
    }
}
