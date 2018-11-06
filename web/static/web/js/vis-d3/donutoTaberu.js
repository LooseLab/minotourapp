"use strict";
// The taxa titles in order to access the results in the AJAX results object
let taxas = ["species", "genus", "family", "order", "classy", "phylum", "superkingdom"];
// the taxa titles we wish to display under the slider
let displayTaxas = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom"];
// Whether drawing for the first time or updating an existing table
let first = true;
// Get the colour scheme, Scale ordinal so colours can be returned the same when given labels
let color = d3.scaleOrdinal(d3.schemeCategory10);
//create data to draw a table
function drawPieTables(countedData, color) {
    //create wrapper function to draw a table
    let tableData = [];
    let obj = {};
    let i;
    // if there is less than 7 members keep the data as one for one table
    // for each data point create data fitting the table
    for (i = 0; i < countedData.length; i++) {
        obj = {
            "Species": countedData[i]["label"],
            "# Matches": countedData[i]["value"],
            "Rank": i + 1,
            "Key": color(countedData[i]["label"])
        };
        tableData.push(obj);
    }
    // if this is the first time drawing the tables for this page
    drawTables("rank-table-container", tableData);
    first = false;
}

// draw the tables using d3
// TODO massively WET code refactor into neater precise code
function drawTables(selection, dataToDraw) {
    // draw the tables using d3
    // get table height
    let columns = ["Species", "# Matches", "Rank", "Key"];
    let Theight = 35;
    let table;
    let thead;
    let tbody;
    let rows;
    let cells;
    // select the table out of one of two sides
    if (first) {
        table = d3.select("." + selection)
            .style("width", "100%").style("height", Theight+"vh").style("overflow-y", "scroll").append("table").attr("class", "table table-hover taxons").style("table-layout", "fixed").style("overflow", "scroll").style("margin-left", "8px");
        // append table head, using the data in the columns array
        thead = table.append('thead').append('tr')
            .selectAll('th')
            .data(columns).enter()
            .append('th')
            .text(function (column) {
                return column;
            });
        // append the table body to the table
        tbody = table.append('tbody').attr("class", selection + "Body");
        rows = tbody.selectAll('tr')
            .data(dataToDraw)
            .enter()
            .append('tr');
        // create cells for the row, one cell for each column
        cells = rows.selectAll('td')
            .data(function (row) {
                return columns.map(function (column) {
                    // if (column !== "Key"){
                    return {column: column, value: row[column]};
                    // }
                });
            })
            .enter()
            .append('td')
            // give the cell the class of the column
            .attr("class", function (d, i) {
                return columns[i];
            })
            .style("background-color", function (d, i) {
                // if the cell contains the key, set the background colour to the rgb value in the data
                if (typeof d.value === "string") {
                    if (d.value.includes("#")) {
                        // d3.select(columns[i]).style()
                        return d.value;
                    }
                    // else set the backgorund to white
                    else {
                        return "white";
                    }
                }
            })
            // fill the cell with the string by setting the inner html
            .html(function (d) {
                if (typeof d.value === "string") {
                    if (!d.value.includes("#")) {
                        return d.value;
                    }
                } else {
                    return d.value;
                }
            });
    } else {
        table = d3.select("." + selection).select("table");
        tbody = table.select("tbody");
        rows = tbody.selectAll("tr");
        // rows.data(dataToDraw).exit().remove();
        rows.data(dataToDraw).enter().append("tr");
        rows.data(dataToDraw).exit().remove();
        rows = tbody.selectAll("tr");
        cells = rows.selectAll('td');
        cells.data(function (row) {
            return columns.map(function (column) {
                // if (column !== "Key"){
                return {column: column, value: row[column]};
                // }
            });
        }).enter().append("td");
        cells = rows.selectAll('td');
        cells.data(function (row) {
            return columns.map(function (column) {
                // if (column !== "Key"){
                return {column: column, value: row[column]};
                // }
            });
        }).attr("class", function (d, i) {
            return columns[i];
        })
        .style("background-color", function (d, i) {
            // if the cell contains the key, set the background colour to the rgb value in the data
            if (typeof d.value === "string") {
                if (d.value.includes("#")) {
                    // d3.select(columns[i]).style()
                    return d.value;
                }
                // else set the backgorund to white
                else {
                    return "white";
                }
            }
        })
        // fill the cell with the string by setting the inner html
        .html(function (d) {
            if (typeof d.value === "string") {
                if (!d.value.includes("#")) {
                    return d.value;
                }
            } else {
                return d.value;
            }
        });

        // give the cell the class of the column
    }

}

function drawDonutRankTable(flowCellId) {
    var selectedBarcode = get_selected_barcode();
    let barcodes = [];
    let flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');
    if (flowcell_selected_tab_input.value !== "Metagenomics"){
        clearInterval(updateDonutTable);
        return;
    }
    $.get("/donut", {flowcellId: flowCellId, visType: "donut", barcode: selectedBarcode}, result => {
        // if there is no data return and try again when interval is up on $interval
        if (result === undefined) {
            return;
        }
        let dataToDraw = result.result;
        let range = $('.input-range'),
             value = $('.taxa-level');
        // what is the value of the current slider? Kingdom, species etc.
        let currentSelectionSlider = value.text();
        // what index is that in the display taxas so we can get the right value from the internal taxas array and results array
        let index = displayTaxas.indexOf(currentSelectionSlider);
        // get the right taxa string, so we can use it as a key o the results object
        let currentlySelectedTaxa = taxas[index];
        // data1 is the data for the currently selected taxa from the results array
        let data1 = dataToDraw[index][currentlySelectedTaxa];
        // get the number of members in this clade used to determine whether we need a second table
        let dataLength = data1.length;
        // if the range slider is changed call the anonymous function to redraw everything
        range.on("input", function () {
            // the selected number for the slider level (0-6)
            let number = this.value;
            // get the right taxa for the key to the results object
            let current_selected_taxa = taxas[number];
            // get the results data array for the currently selected taxa clade
            let sortedData = dataToDraw[number][current_selected_taxa];
            // set the html below the slider to the right level
            value.html(displayTaxas[number]);
            // datalength - how many members ar ein this clade (1-20)
            dataLength = sortedData.length;
            //draw a new donut
            // draw new tables
            drawPieTables(sortedData, color);
        });


        // draw the tables

        drawPieTables(data1, color);
    });
}