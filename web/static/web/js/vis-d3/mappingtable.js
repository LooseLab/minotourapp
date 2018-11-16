function update_mapping_table(flowcellId) {
    let table;
    let thead;
    let tbody;
    let rows;
    let cells;
    let columns = ["Species", "Tax id", "Num. matches", "Prop. classified (%)",
        "Sum. Unique", "Num. mapped",
        "Mapped prop. total (%)", "Danger reads",
        "Red prop. total (%)",
        "Unique Danger reads"
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
        console.log(result);
        result.sort(compare);
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
            .data(result)
            .enter()
            .append('tr').attr("id", function (d) {
            return d.Species.replace(/ /g, "_");
        });
        rows.data(result).exit().remove();
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
        }).attr("class", function (d, i) {
            return columns[i];
        })
            .style("background-color", function (d, i) {
                // if the cell contains the key, set the background colour to the rgb value in the data
                let variable = d3.select("#" + d3.select(this).node().parentNode.id.replace(/ /g, "_"));
                if (d.column === "Num. mapped" && d.value > 0 && !variable.classed("red-alert")) {
                    variable.classed("yellow-alert", false);
                    variable.classed("orange-alert", true);
                }
                else if (d.column === "Danger reads" && d.value > 0) {
                    variable.classed("orange-alert", false);
                    variable.classed("yellow-alert", false);
                    variable.classed("red-alert", true);
                }
                else if (d.column === "Num. matches" && d.value > 0 && !variable.classed("red-alert") && !variable.classed("orange-alert")) {
                    variable.classed("yellow-alert", true);
                }
                else if (d.column === "Num. matches" && d.value=== 0){
                    variable.attr("class", "");
                }
            })
            // fill the cell with the string by setting the inner html
            .html(function (d) {
                return d.value;
            });
    });

}