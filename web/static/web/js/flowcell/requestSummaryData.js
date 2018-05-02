function requestSummaryData (id) {
    /*
     * Request summary by barcode data
     */
    var url = "/api/v1/flowcells/" + id + "/summarybarcode";

    var selectedBarcode = this.getSelectedBarcode();

    $.get(url, (function (data) {

        if (data.length > 0) {

            var basecalled_summary_table = document.querySelector('#basecalled-summary');

            basecalled_summary_table.innerHTML = "";

            data.forEach(function(row) {

                var tr = '<tr>';

                row.forEach(function(column){

                    tr += '<td>' + column + '</td>';
                });

                basecalled_summary_table.innerHTML += tr;
            });
        }

    }).bind(this));
};