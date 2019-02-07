function loadTasksForm() {
    /*
     * Loads create task's form and add event listeners
     */
    // Select the div for the dropdown
    var job_type_select = document.querySelector("#id_job_type");

    var reference_select = document.querySelector("#id_reference");
    // add an event listener for change in the dropdown
    job_type_select.onchange = function (event) {
        // If there is not a value set
        if (!event.target.value) {
            alert('Please select a task.');
        } else {

            if (event.target.value < 0)

                return;

            var selected_job_type;
            // Select the label
            var label = document.querySelector("#select-label");
            // Go through the job types, comapring to see if they are the same as the selected element
            job_type_select.minotour_job_type_list.forEach(function (element) {
                // If they are the same
                if (element["id"] == event.target.value) {
                    //  Create a variable with the selected job type
                    selected_job_type = element;
                    // If it's metagenomics chosen
                    if (element["name"] == "Metagenomics") {
                        // Remove all the reference select options
                        while (reference_select.length > 0) {
                            reference_select.remove(0);
                        }
                        // Fetch the target sets from the database
                        $.get("/api/v1/metagenomics/targetsets", (result, statusText, xhr) => {
                            if (xhr.status !== 204) {
                                label.innerHTML = "Target Sets";
                                for (let i = 0; i < result.length; i++) {
                                    var option = document.createElement("option");
                                    option.value = i;
                                    option.text = result[i];
                                    reference_select.appendChild(option);
                                }
                            }
                        });
                    } else {
                        label.innerHTML = "Reference";
                        while (reference_select.length > 0) {
                            reference_select.remove(0);
                        }
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
                    }
                }

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
