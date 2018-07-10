function requestMinknowMessages(id, data) {

    var url = '/api/v1/flowcells/' + id + '/minknow-messages/?start_time=' + data.runs[0].start_time;

    $.get(url, function (data) {

        stringtowrite = "<div class='table-responsive'><table class='table'><tr><th>Message</th><th>Time</th></tr>";
        for (var i = 0; i < data.length; i++) {
            //stringtowrite=stringtowrite+'<div class="alert alert-info" role="alert">'+data[i].minKNOW_message + ' <p>(<i>' + new Date(data[i].minKNOW_message_timestamp) + '</i>) '+'</div>'
            stringtowrite = stringtowrite + '<tr><td>' + data[i].message + ' </td><td><i>' + new Date(data[i].timestamp) + '</i></td> ' + '</tr>'
        }
        stringtowrite = stringtowrite + '</table></div>';
        document.getElementById('Messages').innerHTML = stringtowrite;
    })
}
