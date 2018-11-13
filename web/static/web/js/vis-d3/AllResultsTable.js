function getTotalReadsTable(flowCellId) {
    var selectedBarcode = get_selected_barcode();
    // Get ttoal reads table updates the total reads table at te bottom of the page
    // Get data from the api
    let flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');
    if (flowcell_selected_tab_input.value !== "Metagenomics") {
        return;
    }
        // Jquery selector
        let table = $(".tableLand");
        // If the table already exists, use the DataTable APi to update in place
        if ($.fn.DataTable.isDataTable(table)) {
            table.DataTable().clear();
            table.DataTable().draw();
        } else {
            // else the databale must be initialised
            table.DataTable({
                    "processing": true,
                    "serverSide": true,
                    "ajax": {
                        "url": '/api/v1/table/',
                        "data": {
                            "flowcell_id": flowCellId
                        }
                    },
                    "columns": [
                        {"data": "barcode_name"},
                        {"data": "superkingdom"},
                        {"data": "phylum"},
                        {"data": "classy"},
                        {"data": "order"},
                        {"data": "family"},
                        {"data": "genus"},
                        {"data": "species"},
                        {"data": "num_matches"},
                        {"data": "proportion_of_classified"},
                    ]
                }
            );
        }
    // updateResultsTable = setTimeout(getTotalReadsTable, 60000, flowCellId, selectedBarcode);
}