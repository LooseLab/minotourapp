function loadTasksForm() {
    /*
     * Loads create task's form and add event listeners
     */

    var job_type_select = document.querySelector("#id_job_type");

    job_type_select.onchange = function (event) {

        if (!event.target.value) {

            alert('Please select a task.');

        } else {

            if (event.target.value < 0)

                return;

            var selected_job_type;

            job_type_select.minotour_job_type_list.forEach(function (element) {

                if (element["id"] == event.target.value)

                    selected_job_type = element;

            });

        }
    };

    var url_task_type_list = "/api/v1/tasktypes/";

    $.get(url_task_type_list, function (dataObj) {

        var data = dataObj['data'];

        while (job_type_select.length > 0) {

            job_type_select.remove(0);
        }

        var option = document.createElement("option");
        option.value = "";
        option.text = "-- select an option --";

        job_type_select.appendChild(option);

        job_type_select.minotour_job_type_list = data; // OO rocks

        for (var i = 0; i < data.length; i++) {

            var option = document.createElement("option");
            option.value = data[i].id;
            option.text = data[i].description;

            job_type_select.appendChild(option);
        }
    });

    var reference_select = document.querySelector("#id_reference");

    var url_reference_list = "/api/v1/reference/";

    $.get(url_reference_list, function (data) {

        let option = document.createElement("option");
        option.value = "";
        option.text = "-- select an option --";

        reference_select.appendChild(option);

        for (var i = 0; i < data.length; i++) {

            option = document.createElement("option");
            option.value = data[i].id;
            option.text = data[i].name;

            reference_select.appendChild(option);
        }
    });


    // $("#post-form-task-create").on("submit", function (event) {
    //     event.preventDefault();
    //
    //
    //     var csrftoken = getCookie('csrftoken');
    //     var csrftoken2 = Cookies.get('csrftoken');
    //
    //     $.ajaxSetup({
    //         beforeSend: function (xhr, settings) {
    //             if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
    //                 xhr.setRequestHeader("X-CSRFToken", csrftoken);
    //             }
    //         }
    //     });
    //
    //     $.ajax({
    //         url: "/api/v1/tasks/",
    //         type: "post",
    //         data: $("#post-form-task-create").serializeArray(),
    //         success: function (json) {
    //             let taskTable = $(".tasktable");
    //             let messages = $("#task_form_messages");
    //             function empty(messages){
    //                 messages.empty();
    //             }
    //             messages.empty();
    //             messages.prepend("<div class=\"alert alert-success\" role=\"alert\">" + json.message + "</div>");
    //             setTimeout(empty, 1500, messages);
    //             taskTable.DataTable().ajax.reload();
    //         },
    //         error: function (json) {
    //
    //             var fields = Object.keys(json['error_messages']);
    //
    //             fields.forEach(function (element) {
    //                 for (var i = 0; i < json['error_messages'][element].length; i++) {
    //                     $("#task_form_messages").prepend("<div class=\"alert alert-danger\" role=\"alert\">" + json['error_messages'][element][i] + "</div>");
    //                 }
    //             });
    //
    //         }
    //     });
    // });
}

function flowcellTaskHistoryTable(flowcellId) {
    // Draw the Task table on the tasks tab showing previous and current analyses
    let table = $(".tasktable");
    // if it's already loaded don't reiniiialise
    if ($.fn.DataTable.isDataTable(table)) {
        table.DataTable().ajax.reload();
    } else {
        table.DataTable({
            ajax: {
                url: '/api/v1/tasks/?search_criteria=flowcell&search_value=' + flowcellId.toString(),
                method: "GET",
            },
            columns: [
                {"data": "id"},
                {"data": "task_type_name"},
                {"data": "last_read"},
                {"data": "read_count"},
                {"data": "running"},
                {"data": "complete"}
            ]
        });
    }
    setTimeout(flowcellTaskHistoryTable, 30000, flowcellId);
}
