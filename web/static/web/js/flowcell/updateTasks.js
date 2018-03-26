function updateTasks(id) {

    var taskstring = "";

    for (var i = 0; i < this.tasks.length; i++) {

        if (this.tasks[i].hasOwnProperty("job_details")) {
            colour = 'bg-green';
            message = 'Reads processed:' + this.tasks[i]["job_details"]["read_count"];
            percentage = 50;
            message2 = 'X% of uploaded reads are processed';
            icon = 'fa fa-refresh fa-spin fa-fw';
            running = true;

        } else {
            colour = 'bg-light-blue';
            message = 'Task Not Running.';
            percentage = 0;
            message2 = "Click to start a " + this.tasks[i]["description"] + ' task.';
            icon = 'fa fa-refresh fa-fw';
            running = false;
        }

        taskstring = this.drawtaskbutton(taskstring, colour, icon, this.tasks[i]["description"], message, percentage, message2, i, this.tasks[i]["long_description"], this.tasks[i]["reference"], this.tasks[i]["name"], this.tasks[i]["transcriptome"]);
    }

    document.getElementById('tasks').innerHTML = taskstring;

    for (var i = 0; i < this.tasks.length; i++) {
        var buttonname = "#button" + this.tasks[i]["name"];
        //console.log("Button name to look for:" + buttonname);
        $(buttonname).click(function (e) {
            var idClicked = e.target.id;
            //console.log(idClicked);
            this.startTask(idClicked.substring(6), id);
        });
    }

    this.liveUpdateTasks(id);

}