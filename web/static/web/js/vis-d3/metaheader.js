"use strict";
// If the table doesn't already exist
let firsty = true;

function metaHeader(flowcellId) {

    // AJAX request for the metadata
    $.get("/centrifuge_metadata", {flowcellId}, result => {

            // Metaheader draws or updates the metadata header at the top of the visualisation page
            // initialise variables
            let table, head, row;
            let target_div, data;
            let flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');

            data = result.result;
            if (flowcell_selected_tab_input.value !== "nav-metagenomics") {
                return;
            }
            if (result.validation){
                console.log("validation");
                target_div = ".meta_taberu";
                d3.select(".meta_taberu_alt").remove();
            } else {
                console.log("no validation");
                target_div = ".meta_taberu_alt";
                d3.select(".meta_taberu").remove();
            }

            if (firsty) {
                // append the table to the page
                console.log(target_div);
                table = d3.select(target_div).append("table").attr("class", "table table-striped").attr("table-layout", "fixed");
                head = table.append("thead");
                row = head.append("tr");
            } else {
                // select the already existing row
                row = d3.select(target_div).select("table").select("tr");
            }
            // If there is no data, return
            if (data[1].value === 0) {
                return;
            }
            // Add any DOM elements for data in the enter selection, which is if there is no th cells already
            row.selectAll("th").data(data).enter().append("th");
            // Update the innerHTML value of the cells
            row.selectAll("th").data(data).html(function (d) {
                return d.key + d.value;
            });
            firsty = false;
        }
    ).fail(function () {
        return;
    });
}