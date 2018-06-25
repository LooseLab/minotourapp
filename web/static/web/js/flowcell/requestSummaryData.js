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

                //console.log(row);

                if (row[0] != "No barcode") {

                    var tr = '<tr>';

                    row.forEach(function (column, index) {

                        if (index == 3 || index == 4 || index == 5 || index == 6) {
                            var text_align = "right";
                        } else {
                            var text_align = "left";
                        }

                        if (index == 6) {
                            column = Math.round(column);
                        }

                        tr += '<td align="' + text_align + '">' + column + '</td>';
                    });

                    basecalled_summary_table.innerHTML += tr;
                }
            });
        }
    }).bind(this));
};