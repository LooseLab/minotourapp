function requestMinknowMessages(id, data) {

    var url = '/api/v1/flowcells/' + id + '/minknow-messages/?start_time=' + data.runs[0].start_time;

    $.get(url, function (data) {

        if (data && data.length > 0) {

            stringtowrite = "";
            stringtowrite = stringtowrite + "<div class='table-responsive'><table class='table'><tr><th>Message</th><th>Time</th></tr>";

            for (var i = 0; i < data.length; i++) {
                stringtowrite = stringtowrite + '<tr><td>' + data[i].message + '</td><td><i>' + new Date(data[i].timestamp) + '</i></td></tr>'
            }

            stringtowrite = stringtowrite + '</table></div>';
            document.getElementById('minknow-messages-div').innerHTML = stringtowrite;

        } else {


            var message = "No information available.";
            var div = document.querySelector("#minknow-messages-div");
            div.innerHTML = message;

        }

    })

}
