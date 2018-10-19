class MTask {

    constructor(flowcell_id, job_type_id, reference_id) {

        this._flowcell_id = flowcell_id;
        this._job_type_id = job_type_id;
        this._reference_id = reference_id;
    }

    get flowcell_id() {

        return this._flowcell_id;
    }

    get job_type_id() {

        return this._job_type_id;
    }

    get reference_id() {

        return this._reference_id;
    }
}