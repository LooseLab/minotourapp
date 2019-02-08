class MTaskController {

    constructor(flowcell_id) {

        this._flowcell_id = flowcell_id;
        this._select_job_type = document.querySelector('#id_job_type');
        this._select_reference = document.querySelector('#id_reference');
        this._task_list = new MTaskList();

        this._message = new MMessage();
        this._messageView = new MMessageView(document.querySelector('#messageView'));
        this._messageView.update(this._message);
    }

    add(event) {

        event.preventDefault();

        let task_new = new MTask(this._flowcell_id, this._select_job_type.value, this._select_reference.value);

        let csrftoken = getCookie('csrftoken');

        let xhr = new XMLHttpRequest();

        xhr.open('POST', '/api/v1/tasks/');

        xhr.setRequestHeader('X-CSRFToken', csrftoken);
        xhr.setRequestHeader("Content-Type", "application/json");
        xhr.setRequestHeader("Accept", "application/json");

        let self = this;

        xhr.onreadystatechange = function() {

            if(this.readyState === XMLHttpRequest.DONE) {

                if(this.status === 200) {

                    self._message.texto = 'Task successfully created!';
                    self._messageView.update(self._message);

                    let taskTable = $(".tasktable");
                    taskTable.DataTable().ajax.reload();

                } else {

                    self._message.texto = 'Something went wrong. Please check the following message. ' + this.responseText;
                    self._messageView.update(self._message);

                }
            }
        };

        xhr.send(JSON.stringify({

            flowcell: task_new.flowcell_id,
            reference: task_new.reference_id,
            job_type: task_new.job_type_id
        }));
    }
}