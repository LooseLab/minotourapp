function requestMinknowMessages(id, data) {

    var url = '/flowcells/' + id + '/minknow_messages_html/?start_time=' + data.runs[0].start_time;

    var run_summaries_table = document.querySelector('#minknow-messages-div');

    $.ajax({
        url: url,
        dataType: "html",
        success: function(data) {
            run_summaries_table.innerHTML = data;
        },
        error: function(e)
        {
            console.log('Error: ' + e);
        }
    });
}
