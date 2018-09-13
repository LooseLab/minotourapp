let updateResultsTable;
function topGetTotalReadsTable(flowCellId) {
    getTotalReadsTable(flowCellId);
    updateResultsTable = setInterval(getTotalReadsTable, 60000, flowCellId);
}

function getTotalReadsTable(flowCellId) {
    // Get ttoal reads table updates the total reads table at te bottom of the page
    // Get data from the api
    let flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');
    if(flowcell_selected_tab_input.value !== "Metagenomics"){
        clearInterval(updateResultsTable);
        console.log("Cleared results table interval");
        return;
    }
    $.get("/table", {flowcellId: flowCellId, visType: "table"}, result => {
        // Jquery selector
        let table = $(".tableLand");
        // If the results of the api query is undefined, there are no results in the DB, so retrun here
        if(result === undefined){
            return;
        }
        // If the table already exists, use the DataTable APi to update in place
        if ($.fn.DataTable.isDataTable(table)) {
            table.DataTable().clear();
            table.DataTable().rows.add(result);
            table.DataTable().draw();
        } else {
            // else the databale must be initialised
            table.DataTable({
                    data: result,
                    columns: [
                        {"data": "superkingdom"},
                        {"data": "phylum"},
                        {"data": "classy"},
                        {"data": "order"},
                        {"data": "family"},
                        {"data": "genus"},
                        {"data": "species"},
                        {"data": "num_matches"},
                        {"data": "prop_species"},
                    ]
                }
            );
        }
    });

}