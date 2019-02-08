"use strict";
// If the table doesn't already exist
let firsty = true;
function metaHeader(flowcellId){
    // Metaheader draws or updates the metadata header at the top of the visualisation page
    // initialise variables
    let table, head, row;
    let flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');
    if(flowcell_selected_tab_input.value !== "nav-metagenomics"){
        return;
    }
    if(firsty){
        // append the table to the page
        table = d3.select(".meta_taberu").append("table").attr("class", "table table-striped").attr("table-layout", "fixed");
        head = table.append("thead");
        row = head.append("tr");
    } else{
        // select the already existing row
        row = d3.select(".meta_taberu").select("table").select("tr");
    }
    // AJAX request for the metadata
    $.get("/centrifuge_metadata", {flowcellId}, result => {
        // If there is no data, return
        if(result[1].value === 0){
            return;
        }
        // Add any DOM elements for data in the enter selection, which is if there is no th cells already
        row.selectAll("th").data(result).enter().append("th");
        // Update the innerHTML value of the cells
        row.selectAll("th").data(result).html(function(d){
            return d.key + d.value;
        });
        }
    ).fail(function(){
        return;
    });
    firsty = false;
}