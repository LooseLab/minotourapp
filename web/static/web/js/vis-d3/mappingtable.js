function update_mapping_table(flowcellId){
    let table;
    let thead;
    let tbody;
    let rows;
    let cells;
    let columns = ["species", "tax_id", "num_matches", "sum_unique","red_reads"];
    let flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');
    if(flowcell_selected_tab_input.value !== "Metagenomics"){
        console.log("Cleared header interval");
        return;
    }
    $.get("/mapped_targets", {flowcellId}, result => {
        console.log(result);

    if(d3.select(".alert-table").classed("has-tabley?")){
        console.log("classy");
        table = d3.select(".alert-table").select("table");
        thead = table.select("thead");
        tbody = table.select("tbody");
    }
    else {
        console.log("not classy");
        table = d3.select(".alert-table").classed("has-tabley?", true).style("width", "100%").append("table").attr("class", "table table-hover");
        thead = table.append('thead').append('tr');
        tbody = table.append('tbody').attr("class", "alert-tbody");

    }

    thead.selectAll('th')
    .data(columns).enter()
    .append('th')
    .text(function (column) {
        console.log(column);
        return column;
    });
    rows = tbody.selectAll('tr');
    rows
        .data(result)
        .enter()
        .append('tr').attr("class", function(d){
            return d.species.replace(/ /g, "_");
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
        if (d.column === "red_reads" && d.value > 0){
            // d3.select("." + )
            d3.select("." + d3.select(this).node().parentNode.className.replace(/ /g, "_")).style("background-color", "red");
        } else {
            d3.select("." + d3.select(this).node().parentNode.className.replace(/ /g, "_")).style("background-color", "yellow");
        }
    })
    // fill the cell with the string by setting the inner html
    .html(function (d) {
        return d.value;
    });



    console.log(flowcellId);
    });

}