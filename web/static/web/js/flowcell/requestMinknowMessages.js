function requestMinknowMessages(id) {

    console.log("Requesting Minknow messages...");
    
    if (!this.rundetails) {
        return;
    }

    var url_sincemessages = this.rundetails[0]["minION"] + 'messagessince/' + this.rundetails[0]["minKNOW_start_time"] + '/' + this.lasttime.toISOString() + "/";

    console.log(url_sincemessages);

    $.get(url_sincemessages, function (data) {
        console.log(data);
        stringtowrite = '<table class="table table-condensed"><tr><th>Message</th><th>Time</th></tr>';
        for (var i = 0; i < data.length; i++) {
            //stringtowrite=stringtowrite+'<div class="alert alert-info" role="alert">'+data[i].minKNOW_message + ' <p>(<i>' + new Date(data[i].minKNOW_message_timestamp) + '</i>) '+'</div>'
            stringtowrite = stringtowrite + '<tr><td>' + data[i].minKNOW_message + ' </td><td><i>' + new Date(data[i].minKNOW_message_timestamp) + '</i></td> ' + '</tr>'
        }
        stringtowrite = stringtowrite + '</table>';
        document.getElementById('Messages').innerHTML = stringtowrite;
    })
}
