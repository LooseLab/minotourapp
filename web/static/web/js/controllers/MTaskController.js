class MTaskController {

    constructor(flowcell_id) {

        this._flowcell_id = flowcell_id;
        this._select_job_type = document.querySelector('#id_job_type');
        this._select_reference = document.querySelector('#id_reference');
        this._task_list = new MTaskList();
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

        xhr.onreadystatechange = function() {

            if(this.readyState === XMLHttpRequest.DONE) {

                console.log('Response status: ' +  this.status);
                console.log('Response text: ' + this.responseText);

                if(this.status === 200) {

                    let taskTable = $(".tasktable");

                    let messages = $("#task_form_messages");

                    function empty(messages){
                        messages.empty();
                    }

                    messages.empty();

                    messages.prepend("<div class=\"alert alert-success\" role=\"alert\">Success</div>");

                    setTimeout(empty, 1500, messages);

                    taskTable.DataTable().ajax.reload();
                }
            }
        };

        console.log(JSON.stringify({

            flowcell: task_new.flowcell_id,
            reference: task_new.reference_id,
            job_type: task_new.job_type_id
        }));

        xhr.send(JSON.stringify({

            flowcell: task_new.flowcell_id,
            reference: task_new.reference_id,
            job_type: task_new.job_type_id
        }));
    }
}