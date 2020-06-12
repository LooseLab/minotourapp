function loadTasksForm() {
    /** TODO rewrite
     * @function
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
                        document.getElementById("id_reference").disabled = false;

                        while (reference_select.options.length > 0) {
                            reference_select.remove(0);
                        }
                        // Fetch the target sets from the database
                        $.get("/api/v1/metagenomics/targetsets", (result, statusText, xhr) => {
                            if (xhr.status !== 204) {
                                // Change the second dropdown label text to target sets
                                label.innerHTML = "Target Sets";
                                // Create an option element for the placeholder for the dropdown
                                let option = document.createElement("option");

                                option.value = "";
                                // set placeholder text
                                option.text = "-- select an option --";
                                // add the placeholder option to the second drop down
                                reference_select.appendChild(option);
                                // for all the target sets
                                for (let i = 0; i < result.length; i++) {
                                    // create options for each available target set
                                    option = document.createElement("option");
                                    option.value = i;
                                    option.text = result[i];
                                    reference_select.appendChild(option);
                                }
                            }
                        });
                        // if it's a minimap2 task, make the reference select box to be usable.
                    } else if (element["name"] === "Minimap2"){
                        document.getElementById("id_reference").disabled = false;
                    } else {
                        // if it's not for the metagenomics task we label the second dropdown reference
                        document.getElementById("id_reference").disabled = true;
                        label.innerHTML = "Reference";
                        // remove all the current references
                        while (reference_select.options.length > 0) {
                            reference_select.remove(0);
                        }
                        // Url to fetch available references
                        var url_reference_list = "/api/v1/reference/";

                        $.get(url_reference_list, function (data) {
                            // Create a placeholder option
                            let option = document.createElement("option");

                            option.value = "";
                            // set the placeholder option text to
                            option.text = "-- select an option --";
                            // add the placeholder option to the select element
                            reference_select.appendChild(option);
                            // for each available reference
                            for (var i = 0; i < data.length; i++) {
                                // create an option element
                                option = document.createElement("option");
                                // set the value to the reference id
                                option.value = data[i].id;
                                // set the value text to the reference name
                                option.text = data[i].name;
                                // add the child option
                                reference_select.appendChild(option);
                            }
                        });
                    }
                }

            });

        }
    };
    // url to fetch all the available task types
    var url_task_type_list = "/api/v1/tasktypes/";
    // fetch all available tasks types
    $.get(url_task_type_list, function (dataObj, statusText, xhr) {
        if (xhr.status != 200) console.log(`Error ${statusText}`);
        // get the data
        var data = dataObj['data'];
        // create option placeholder
        var option = document.createElement("option");

        let alreadyPresentJobs = [...$("#id_job_type")[0].children];

        if (!alreadyPresentJobs.length){
            option.value = "";
            // set placeholder text
            option.text = "-- select an option --";
            // add the option to the job select element
            job_type_select.appendChild(option);
            // set this attribute on the elect box
            job_type_select.minotour_job_type_list = data; // OO rocks
        }

        // for each of the elements
        alreadyPresentJobs.forEach(function (presentJob) {
            // if the chromsome in already present chromsome is in the list of chromosomes fetched from the server
            if (data.some(e => e.description === presentJob.textContent)){
                // get the index,
                let index = data.map(function(x) {return x.description; }).indexOf(presentJob.textContent);

                // remove it from the freshly fetched chromosomes
                data.splice(index, 1);
            }
        });


        if (!data.length) {return;}

        // for each job type
        for (var i = 0; i < data.length; i++) {
            // create an option
            var option = document.createElement("option");
            // Set th option value as the job type id
            option.value = data[i].id;
            // Set the text to the description
            option.text = data[i].description;
            // Add the job type option the the select element
            job_type_select.appendChild(option);
        }
    });
    // Do the same with the reference

    if (document.getElementById("id_job_type").value !== "4"){
        document.getElementById("id_reference").disabled=true;
    }

    var url_reference_list = "/api/v1/reference/";

    $.get(url_reference_list, data => {

        let option = document.createElement("option");

        let references = data;

        let alreadyPresentRefs = [...$("#id_reference")[0].children];

        if (!alreadyPresentRefs.length){
            option.value = "";
            option.text = "-- select an option --";
            reference_select.appendChild(option);
        }
        // for each of the elements
        alreadyPresentRefs.forEach(function (presentRef) {
            // if the chromsome in already present chromsome is in the list of chromosomes fetched from the server
            if (references.some(e => e.name === presentRef.textContent)) {
                // get the index,
                let index = references.map(function(x) {return x.name; }).indexOf(presentRef.textContent);
                // remove it from the freshly fetched chromosomes
                references.splice(index, 1);
            }
        });

        if (!references.length) {return;}

        for (var i = 0; i < references.length; i++) {

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
                {"data": "read_count"},
                {"data": "running"},
                {"data": "complete"},
                {"data": "reference_name"},
                {"data": "server_initiated"},
                {
                    "data": null,
                    "orderable": false,
                    "width": "15%",
                    "render": (data, type, full) => {
                        if (!data["server_initiated"]) {
                            return [`<a class="btn icon-task"  id="pause_${data.id}" onclick="mTaskController.performActionOnTask(event, ${data.id}, 2)" ><i class="fa fa-${data.icon} task-icon"></i> ${data.iconText} </a>·
                                <a class="btn icon-task" id="reset_${data.id}" onclick="mTaskController.performActionOnTask(event, ${data.id}, 1)"><i class="fa fa-recycle task-icon"></i> Reset </a>·
                                <a class="btn icon-task" id="delete_${data.id}" onclick="mTaskController.performActionOnTask(event, ${data.id}, 3)"><i class="fa fa-times task-icon"></i> Delete </a>`,];
                        }else{
                            return [`<a></a>`]
                        }
                    }
                }
            ]
        });
    }
}
