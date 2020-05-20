function requestMinknowMessages(id, data) {
    // Request minKnow messages for the live events page and flowcell summary page
    var url = '/api/v1/flowcells/' + id + '/minknow-messages-html/';
    // Select the message table node in the DOM
    var minknowMessagesDiv = document.querySelector('#minknow-messages-div');
    console.log(data);
    $.ajax({
        url: url,
        dataType: "html",
        // Append the rendered HTML to the DOM node we selected above, django is returning a fully rendered HTML template
        success: function(data) {
            minknowMessagesDiv.className="";
            minknowMessagesDiv.innerHTML = data;
        },
        // Send the get request params, start time of the first run, end time of the final run
        // TODO needs testing with multiple runs, as I'm not sure runs will be in the right order

        // If there is no last read time (end time)
        data:{},
        error: function(e)
        {
            console.log('Error: ' + e);
        }
    });
}
