function requestSummaryData(id) {
    /*
     * Request flowcell summary
     * This response is a HTML code that is appended to the page (no json)
     */
    var url = "/flowcells/" + id + "/summary_html";

    var selected_barcode = get_selected_barcode();

    var basecalled_summary_table = document.querySelector('#basecalled-summary');

    $.ajax({
        url: url,
        dataType: "html",
        success: function(data) {
            basecalled_summary_table.innerHTML = data;
        },
        error: function(e)
        {
            console.log('Error: ' + e);
        }
    });
};
