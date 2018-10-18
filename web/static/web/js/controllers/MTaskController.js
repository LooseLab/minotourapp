class MTaskController {

    constructor(flowcell_id) {

        this._flowcell_id = flowcell_id;
        this._select_job_type = document.querySelector('#id_job_type');
        this._select_reference = document.querySelector('#id_reference');
    }

    add() {

        let task_new = new MTask(this._flowcell_id, this._select_job_type.value, this._select_reference.value);
        console.log(task_new);
    }
}