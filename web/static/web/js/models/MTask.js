class MTask {
    // a class for a Job to be created - stores all the attributes we might need
    constructor(flowcell_id, job_type_id, reference_id, taret_set_id) {
        // construct all the attributes we might need
        this._flowcell_id = flowcell_id;
        this._job_type_id = job_type_id;
        this._reference_id = reference_id;
        this._target_set_id = taret_set_id;
    }
    // getter for flowcell id
    get flowcell_id() {

        return this._flowcell_id;
    }
    // getter for job type id
    get job_type_id() {

        return this._job_type_id;
    }
    // getter for the reference ID
    get reference_id() {

        return this._reference_id;
    }
    // getter for the target_set id
    get target_set_id() {

        return this._target_set_id;
    }
}