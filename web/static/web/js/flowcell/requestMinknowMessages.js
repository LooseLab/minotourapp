function requestMinknowMessages(id, data) {
    // Request minKnow messages for the live events page and flowcell summary page
    var url = '/flowcells/' + id + '/minknow_messages_html/';
    // Select the message table node in the DOM
    var minknowMessagesDiv = document.querySelector('#minknow-messages-div');
    console.log(data);
    $.ajax({
        url: url,
        dataType: "html",
        // Append the rendered HTML to the DOM node we selected above, django is returning a fully rendered HTML template
        success: function(data) {
            minknowMessagesDiv.innerHTML = data;
        },
        // Send the get request params, start time of the first run, end time of the final run
        // TODO needs testing with multiple runs, as I'm not sure runs will be in the right order

        // If there is no last read time (end time)
        data:{start_time: data.runs[0].start_time, end_time:data.runs[data.runs.length-1].last_read_time},
        error: function(e)
        {
            console.log('Error: ' + e);
        }
    });
}
