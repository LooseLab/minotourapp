function requestMinknowMessages(id, data) {
    // Request minKnow messages for the live events page
    var url = '/flowcells/' + id + '/minknow_messages_html/';
    // Select the table node in the DOM
    var run_summaries_table = document.querySelector('#minknow-messages-div');

    $.ajax({
        url: url,
        dataType: "html",
        // Append the rendered HTML to the DOM node we selected above, django is returning a fully rendered HTML template
        success: function(data) {
            run_summaries_table.innerHTML = data;
        },
        // Send the get request params, start time of the first run, end time of the final run
        // TODO needs testing with multiple runs, as I'm not sure runs will be in the right order
        data:{start_time: data.runs[0].start_time, end_time:data.runs[data.runs.length-1].last_read_time},
        error: function(e)
        {
            console.log('Error: ' + e);
        }
    });
}
