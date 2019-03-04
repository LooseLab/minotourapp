function revealMetagenomicsPage(){
    d3.select("#loading-sign").transition().duration(1000).style("opacity", 0);

    setTimeout(function () {
        $("body").addClass("loaded sidebar-collapse");
        d3.select("#loading-sign").style("display", "none");
        d3.select(".vis-container").style("display", "contents");
    }, 1500);
}

function update_mapping_table(flowcellId) {
    let table;
    let thead;
    let tbody;
    let rows;
    let cells;
    let data;
    let row_count = 0;
    let columns = ["Alert level", "Species", "Num. matches", "Prop. classified (%)",
        "Num. mapped",
        "Mapped prop. matches (%)", "Target reads",
        "Target prop. mapped (%)",
         ];
    let flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');
    let barcode = get_selected_barcode();

    // Order the results correctly
    function compare(a, b) {
        if (a.Species < b.Species)
            return -1;
        if (a.Species > b.Species)
            return 1;
        return 0;
    }

    if (flowcell_selected_tab_input.value !== "nav-metagenomics") {
        return;
    }
    $.get("/mapped_targets", {flowcellId, barcode}, (result, statusText, xhr) => {
        if (xhr.status == 204){
            let alertTable = d3.select(".alert-table-complex");
            d3.select(".button-row").remove();
            d3.select("#demo").classed("in", true);
            alertTable.remove();
            revealMetagenomicsPage();
            return;
        }
        if (result.table === undefined){
            let alertTable = d3.select(".alert-table");
            alertTable.selectAll("*").remove();
            alertTable.append("text").text("No data to display");
            alertTable.classed("no-data", true);
            revealMetagenomicsPage();
            return;
        } else {
            d3.select(".alert-table").classed("no-data", false).select("text").remove();
            revealMetagenomicsPage();
        }

        data = result.table;
        data.sort(compare);
        if (d3.select(".alert-table").classed("has-tabley?")) {
            table = d3.select(".alert-table").select("table");
            thead = table.select("thead");
            tbody = table.select("tbody");
        }
        else {
            table = d3.select(".alert-table").classed("has-tabley?", true).style("width", "100%").append("table").attr("class", "table map-alert");
            thead = table.append('thead').append('tr');
            tbody = table.append('tbody').attr("class", "alert-tbody");

        }

        thead.selectAll('th')
            .data(columns).enter()
            .append('th')
            .text(function (column) {
                return column;
            });
        rows = tbody.selectAll('tr');
        rows
            .data(data)
            .enter()
            .append('tr').attr("id", function (d) {
            return d.Species.replace(/ /g, "_");
        });
        rows.data(data).exit().remove();
        rows = tbody.selectAll('tr');
        cells = rows.selectAll('td');
        cells.data(function (row) {
            return columns.map(function (column) {
                // if (column !== "Key"){
                return {column: column, value: row[column]};
                // }
            });
        })
            .enter()
            .append('td');
        cells = rows.selectAll("td");
        cells.data(function (row) {
            return columns.map(function (column) {
                // if (column !== "Key"){
                return {column: column, value: row[column]};
                // }
            });
        }).attr("id", function (d, i) {
            if (d.column === "Alert level"){
                row_count += 1;
            }
            return columns[i].replace(" ", "_") + "_" + row_count.toString();
        })
        .style("background-color", function (d, i) {
            // if the cell contains the key, set the background colour to the rgb value in the data
            let variable = d3.select("#" + d3.select(this).node().parentNode.children[0].id);
            if (d.column === "Num. mapped" && d.value > 0) {
                variable.classed("green-alert", false);
                variable.classed("yellow-alert", false);
                variable.classed("orange-alert", true);
            }
            else if (d.column === "Target reads" && d.value > 0) {
                variable.classed("green-alert", false);
                variable.classed("orange-alert", false);
                variable.classed("yellow-alert", false);
                variable.classed("red-alert", true);
            }
            else if (d.column === "Num. matches" && d.value > 0) {
                variable.classed("green-alert", false);
                variable.classed("yellow-alert", true);
            }
            else if (d.column === "Num. matches" && d.value=== 0){
                variable.attr("class", "green-alert");
            }
        })
        // fill the cell with the string by setting the inner html
        .html(function (d) {
            return d.value;
        });
    });

}