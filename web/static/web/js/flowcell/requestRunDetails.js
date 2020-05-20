function requestRunDetails(id) {

    var url = "/flowcells/" + id + "/run_summaries_html/";

    var run_summaries_table = document.querySelector('#run-summaries-div');

    $.ajax({
        url: url,
        dataType: "html",
        success: function(data) {
            run_summaries_table.className = "";
            run_summaries_table.innerHTML = data;
        },
        error: function(e)
        {
            console.log('Error: ' + e);
        }
    });

    var url = "/flowcells/" + id + "/flowcell_run_basecalled_summary_html/";

    var flowcell_run_basecalled_summary_html_div = document.querySelector('#flowcell_run_basecalled_summary_html_div');

    $.ajax({
        url: url,
        dataType: "html",
        success: function(data) {
            flowcell_run_basecalled_summary_html_div.className = "";
            flowcell_run_basecalled_summary_html_div.innerHTML = data;
        },
        error: function(e)
        {
            console.log('Error: ' + e);
        }
    });


};
