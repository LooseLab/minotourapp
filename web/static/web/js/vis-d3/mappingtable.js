function update_mapping_table(flowcellId) {
    let table;
    let thead;
    let tbody;
    let rows;
    let cells;
    let row_count = 0;
    let columns = ["Alert level", "Species", "Num. matches", "Prop. classified (%)",
        "Sum. Unique", "Num. mapped",
        "Mapped prop. total (%)", "Target reads",
        "Red prop. total (%)",
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
    $.get("/mapped_targets", {flowcellId, barcode}, result => {
        // Set the barcode tabs to the highest alert level in their contents
        let alertLevels = {0: "green-alert-tab", 1 : "yellow-alert-tab", 2 : "orange-alert-tab", 3 : "red-alert-tab"};
        let tab_level;
        for (let i = 0; i < result.tabs.length; i++){
            tab = result.tabs[i];
            console.log(tab);
            tab_level = alertLevels[tab["value"]];
            console.log(tab_level);
            d3.select("."+tab["key"]).classed(tab_level, true);
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