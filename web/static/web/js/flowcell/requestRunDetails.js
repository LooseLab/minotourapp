function requestRunDetails(id) {

    var url = "/flowcells/" + id + "/run_summaries_html/";

    var run_summaries_table = document.querySelector('#run-summaries-div');

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
};
