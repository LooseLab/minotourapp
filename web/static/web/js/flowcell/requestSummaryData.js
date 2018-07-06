function requestSummaryData(id) {
    /*
     * Request summary by barcode data
     */
    var url = "/api/v1/flowcells/" + id + "/summarybarcode";

    var selectedBarcode = this.getSelectedBarcode();

    $.get(url, (function (data) {

        if (data.length > 0) {

            var basecalled_summary_table = document.querySelector('#basecalled-summary');

            basecalled_summary_table.innerHTML = "";

            data.forEach(function (row) {

                var tr = '<tr>';
                tr += '<td>' + row['barcode_name'] + '</td>';
                tr += '<td>' + row['read_type_name'] + '</td>';
                tr += '<td>' + row['status'] + '</td>';
                tr += '<td align="right">' + row['read_count'] + '</td>';
                tr += '<td align="right">' + row['total_length'] + '</td>';
                tr += '<td align="right">' + row['total_length'] + '</td>';
                tr += '<td align="right">' + row['total_length'] + '</td>';

                // if (index == 6) {
                //     column = Math.round(column);
                // }

                basecalled_summary_table.innerHTML += tr;

            });
        }
    }).bind(this));
};