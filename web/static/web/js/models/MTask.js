class MTask {

    constructor(flowcell_id, task_type_id, reference_id) {

        this._flowcell_id = flowcell_id;
        this._task_type_id = task_type_id;
        this._reference_id = reference_id;
    }

    get flowcell_id() {

        return this._flowcell_id;
    }
}