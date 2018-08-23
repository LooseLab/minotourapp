function liveUpdateTasks(id) {

    var url = "/api/v1/flowcells/" + id + "/tasks/";

    $.get(url, (function (tasks) {

        for (var i = 0; i < tasks.length; i++) {

            if (tasks[i].hasOwnProperty("job_details") && this.summary) {

                message = 'Reads processed:' + tasks[i]["job_details"]["read_count"] + "/" + this.summary["All reads"]["Template"]["read_count"]["all"]["data"][0];
                percentage = Math.round((tasks[i]["job_details"]["read_count"] / this.summary["All reads"]["Template"]["read_count"]["all"]["data"][0] * 100) * 100) / 100;
                message2 = percentage + '% of uploaded reads are processed';

                $('#' + self.tasks[i]["name"] + '-message').text(message);
                $('#' + self.tasks[i]["name"] + '-percentage').text(percentage);
                $('div#' + self.tasks[i]["name"] + '-percentage').width(percentage + '%');
                $('#' + self.tasks[i]["name"] + '-message2').text(message2);

            }

        }

    }).bind(this));
};

function startTask(description, reference) {

    var e = document.getElementById(description);

    if (e != null) {

        var strUser = e.options[e.selectedIndex].value;

    } else {

        var strUser = "null";

    }

    $.ajaxSetup({
        beforeSend: function (xhr, settings) {
            if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
                xhr.setRequestHeader("X-CSRFToken", getCookie('csrftoken'));
            }
        }
    });


    var url_mask = "/api/v1/flowcells/12345/settask/".replace(/12345/, reference.toString());

    $.ajax({
        "type": "POST",
        "dataType": "json",
        "url": url_mask,
        "data": {
            "job": description,
            "reference": strUser,
        },
        "beforeSend": function (xhr, settings) {
            //console.log("before send");
            $.ajaxSettings.beforeSend(xhr, settings);
        },
        "success": function (result) {
            //console.log(result);
            $(".modal.in").modal("hide");
            this.requestTasks(reference);

        },
        error: function (XMLHttpRequest, textStatus, errorThrown) {
            //console.log("Status: " + textStatus);
            //alert("Error: " + errorThrown);
        }

    })
};


function drawtaskbutton(taskstring, colour, icon, description, message, percentage, message2, i, long_description, reference, name, transcriptome) {

    taskstring = taskstring + '<div class="col-md-4">';
    taskstring = taskstring + '<button type="button" class="info-box ' + colour + '" data-toggle="modal" data-target="#taskmodal' + i + '">';
    taskstring = taskstring + '<div >';
    taskstring = taskstring + '<span class="info-box-icon"><i class="' + icon + '"></i></span>';
    taskstring = taskstring + '<div class="info-box-content">';
    taskstring = taskstring + '<span class="info-box-text">' + description + '</span>';

    taskstring = taskstring + '<span class="info-box-number" id="' + name + '-message">' + message + '</span>';
    taskstring = taskstring + '<div class="progress">';
    taskstring = taskstring + '<div class="progress-bar" id="' + name + '-percentage" style="width: ' + percentage + '%"></div>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '<span class="progress-description" id="' + name + '-message2" >';
    taskstring = taskstring + message2;
    taskstring = taskstring + '</span>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '</button>';
    taskstring = taskstring + '<div id="taskmodal' + i + '" class="modal fade" role="dialog">';
    taskstring = taskstring + '<div class="modal-dialog">';
    taskstring = taskstring + '<div class="modal-content">';
    taskstring = taskstring + '<div class="modal-header">';
    taskstring = taskstring + '<button type="button" class="close" data-dismiss="modal">&times;</button>';
    taskstring = taskstring + '<h4 class="modal-title">' + description + '</h4>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '<div class="modal-body">';
    taskstring = taskstring + '<p>' + long_description + '</p>';

    if (reference == true) {
        taskstring = taskstring + "<p>Select Reference:</p>"
        taskstring = taskstring + '<select id="' + name + '">';
        //console.log(this.references);
        if (this.references) {
            for (var j = 0; j < this.references.length; j++) {
                if (this.references[j]['transcripts'] == transcriptome) {
                    taskstring = taskstring + '<option>' + this.references[j]["reference_name"] + '</option>';
                }
            }
        }
        taskstring = taskstring + '</select>';
    }

    taskstring = taskstring + '</div>';
    taskstring = taskstring + '<div class="modal-footer">';
    taskstring = taskstring + '<button type="button" class="btn btn-default" id="button' + name + '">Go</button>';
    taskstring = taskstring + '<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '</div>';

    return taskstring;
};

function updateTasks(tasks) {

    var taskstring = "";

    for (var i = 0; i < tasks.length; i++) {

        if (tasks[i].hasOwnProperty("job_details")) {

            colour = 'bg-green';
            message = 'Reads processed:' + tasks[i]["job_details"]["read_count"];
            percentage = 50;
            message2 = 'X% of uploaded reads are processed';
            icon = 'fa fa-refresh fa-spin fa-fw';
            running = true;

        } else {

            colour = 'bg-light-blue';
            message = 'Task Not Running.';
            percentage = 0;
            message2 = "Click to start a " + tasks[i]["description"] + ' task.';
            icon = 'fa fa-refresh fa-fw';
            running = false;

        }

        taskstring = drawtaskbutton(
            taskstring,
            colour,
            icon,
            tasks[i]["description"],
            message,
            percentage,
            message2,
            i,
            tasks[i]["long_description"],
            tasks[i]["reference"],
            tasks[i]["name"],
            tasks[i]["transcriptome"]
        );
    }

    document.getElementById('tasks').innerHTML = taskstring;

    var flowcell_id = get_selected_flowcell();

    for (var i = 0; i < tasks.length; i++) {

        var buttonname = "#button" + tasks[i]["name"];

        $(buttonname).click(function (e) {

            var idClicked = e.target.id;

            startTask(idClicked.substring(6), flowcell_id);
        });
    }

    liveUpdateTasks(flowcell_id);
}

function requestTasks(id) {

    var url = "/api/v1/flowcells/" + id + "/tasks/";

    $.get(url, (function (data) {
        updateTasks(data); //really this only needs to run once!
    }).bind(this));
}

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

            // console.log('You selected task ' + selected_job_type["description"] + '.');
        }
    };

    var url_task_type_list = "/api/v1/tasktypes/";

    $.get(url_task_type_list, (function (dataObj) {

        var data = dataObj['data'];

        while (job_type_select.length > 0) {

            job_type_select.remove(0);
        }

        var option = document.createElement("option");
        option.value = -1;
        option.text = "-- select an option --";

        job_type_select.appendChild(option);

        job_type_select.minotour_job_type_list = data; // OO rocks

        for (var i = 0; i < data.length; i++) {

            var option = document.createElement("option");
            option.value = data[i].id;
            option.text = data[i].description;

            job_type_select.appendChild(option);
        }
    }));

    var reference_select = document.querySelector("#id_reference");

    var url_reference_list = "/api/v1/reference/";

    $.get(url_reference_list, (function (data) {

        var option = document.createElement("option");
        option.value = -1;
        option.text = "-- select an option --";

        reference_select.appendChild(option);

        for (var i = 0; i < data.length; i++) {

            var option = document.createElement("option");
            option.value = data[i].id;
            option.text = data[i].reference_name;

            reference_select.appendChild(option);
        }
    }));

    /*var btn_task_create = document.querySelector("#btn_task_create");

    btn_task_create.addEventListener("click", function(event) {

        event.preventDefault();

        var task_type_select = document.querySelector("#task_type_id");
        var reference_select = document.querySelector("#reference_id");

        console.log('task_type_select: ' + task_type_select.value);

    });*/

    $("#post-form-task-create").on("submit", function (event) {
        event.preventDefault();
        console.log("form submitted");

        console.log($(this).serializeArray());

        var csrftoken = getCookie('csrftoken');
        var csrftoken2 = Cookies.get('csrftoken');

        $.ajaxSetup({
            beforeSend: function (xhr, settings) {
                if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
                    xhr.setRequestHeader("X-CSRFToken", csrftoken);
                }
            }
        });

        $.ajax({
            url: "/api/v1/tasks/",
            type: "post",
            data: $("#post-form-task-create").serializeArray(),
            success: function (json) {
                console.log(json);
                console.log('success');
                $("#task_form_messages").prepend("<div class=\"alert alert-success\" role=\"alert\">" + json.message + "</div>");
            },
            error: function (json) {
                console.log(json);
                console.log('error');

                var fields = Object.keys(json['error_messages']);
                console.log(fields);

                fields.forEach(function(element) {
                    console.log(element);
                    for (var i = 0; i < json['error_messages'][element].length; i++) {
                        $("#task_form_messages").prepend("<div class=\"alert alert-danger\" role=\"alert\">" + json['error_messages'][element][i] + "</div>");
                    }
                });

            }
        });
    });
}