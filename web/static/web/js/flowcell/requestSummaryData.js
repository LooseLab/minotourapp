function requestSummaryData(id) {
    /*
     * Request flowcell summary
     * This response is a HTML code that is appended to the page (no json)
     */
    var url = "/flowcells/" + id + "/summary_html";

    var selectedBarcode = this.getSelectedBarcode();

    var basecalled_summary_table = document.querySelector('#basecalled-summary');

    console.log('requesting summary data');
    $.ajax({
        url: url,
        dataType: "html",
        success: function(data) {
            basecalled_summary_table.innerHTML = data;
            console.log('request returned success');
        },
        error: function(e)
        {
            alert('Error: ' + e);
        }
    });
};