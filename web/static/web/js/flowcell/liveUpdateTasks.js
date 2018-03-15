function liveUpdateTasks(id) {

    var url = "/api/v1/flowcells/" + id + "/tasks/";

    $.get(url, (function (tasks) {

        for (var i = 0; i < tasks.length; i++) {

            if (tasks[i].hasOwnProperty("job_details")) {

                console.log('> summary');
                console.log(this.summary);
                console.log('< summary')

                message = 'Reads processed:' + tasks[i]["job_details"]["read_count"] + "/" + this.summary["All reads"]["Template"]["read_count"]["all"]["data"][0];
                percentage = Math.round((tasks[i]["job_details"]["read_count"] / this.summary["All reads"]["Template"]["read_count"]["all"]["data"][0] * 100) * 100) / 100;
                message2 = percentage + '% of uploaded reads are processed';
                //console.log(message);
                $('#' + self.tasks[i]["name"] + '-message').text(message);
                $('#' + self.tasks[i]["name"] + '-percentage').text(percentage);
                $('div#' + self.tasks[i]["name"] + '-percentage').width(percentage + '%');
                $('#' + self.tasks[i]["name"] + '-message2').text(message2);

            }

        }

    }).bind(this));

};