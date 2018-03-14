function liveUpdateTasks(id) {
    var url = "/api/v1/flowcells/" + id + "/tasks/";
    $.get(url, function (data) {
        var tasks = [];
        for (var i = 0; i < data.length; i++) {
            tasks.push(data[i]);
        }
        self.tasks = tasks;
        for (var i = 0; i < self.tasks.length; i++) {
            if (self.tasks[i].hasOwnProperty("job_details")) {
                message = 'Reads processed:' + self.tasks[i]["job_details"]["read_count"] + "/" + self.summary["All reads"]["Template"]["read_count"]["all"]["data"][0];
                percentage = Math.round((self.tasks[i]["job_details"]["read_count"] / self.summary["All reads"]["Template"]["read_count"]["all"]["data"][0] * 100) * 100) / 100;
                message2 = percentage + '% of uploaded reads are processed';
                //console.log(message);
                $('#' + self.tasks[i]["name"] + '-message').text(message);
                $('#' + self.tasks[i]["name"] + '-percentage').text(percentage);
                $('div#' + self.tasks[i]["name"] + '-percentage').width(percentage + '%');
                $('#' + self.tasks[i]["name"] + '-message2').text(message2);
            }
        }
    })
};