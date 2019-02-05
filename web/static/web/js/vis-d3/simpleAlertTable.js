function compare(a, b) {
    if (a["Potential threats"] < b["Potential threats"])
        return -1;
    if (a["Potential threats"] > b["Potential threats"])
        return 1;
    return 0;
}

function draw_simple_table(flowcellId) {
    let table;
    let thead;
    let tbody;
    let rows;
    let cells;
    let columns = ["1 in 5000000", "1 in 2000000", "1 in 1000000", "1 in 100000", "1 in 50000", "1 in 10000",
        "1 in 1000", "Potential threats", "Detected"];
    let data;
    let limit;
    let flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');
    let barcode = get_selected_barcode();
    if (flowcell_selected_tab_input.value !== "nav-metagenomics") {
        return;
    }
    $.get("/api/v1/metagenomics/alerts", {flowcellId, barcode}, (result, statusText, xhr) => {
        if (xhr.status == 204){
            d3.select(".alert-table-simple").remove();
            return;
        }
        data = result.table;
        limit = result.conf_limit;
        data.sort(compare);
        if (d3.select(".simple-alert-table").classed("has-tabley?")) {
            table = d3.select(".simple-alert-table").select("table");
            thead = table.select("thead");
            tbody = table.select("tbody");
        }
        else {
            table = d3.select(".simple-alert-table").classed("has-tabley?", true)
                .style("width", "100%").append("table").attr("class", "table map-alert simple")
                .style("table-layout", "fixed").style("border-bottom", "1px solid #f4f4f4");
            thead = table.append('thead').append('tr');
            tbody = table.append('tbody').attr("class", "alert-tbody");
        }
        thead.selectAll('th')
            .data(columns).enter()
            .append('th').style("text-align", "center")
            .attr("id", function (column) {
                return column.replace(" ", "_");
            }).attr("class", function (d, i) {
            return i;
            }).text(function (column) {
                return column;
            });

        d3.select("#Not_detected").html("Not detected at 1 in " + limit);
        rows = tbody.selectAll('tr');
        rows
            .data(data)
            .enter()
            .append('tr').attr("id", function (d) {
            return d["Potential threats"].replace(/ /g, "_");
        });
        rows.data(data).exit().remove();
        rows = tbody.selectAll('tr');
        cells = rows.selectAll('td');
        cells.data(function (row) {
            return columns.map(function (column) {
                return {column: column, value: row[column]};
            });
        })
            .enter()
            .append('td');
        cells = rows.selectAll("td");
        cells.data(function (row) {
            console.log(row);
            return columns.map(function (column) {
                return {
                    column, value: row[column], read_count: row.read_count,
                    detected: row.Detected, species: row["Potential threats"], index: row.column_index,
                    conf_limit: row.conf_limit
                };
            });
        }).attr("id", function (d, i) {
            return columns[i].replace(" ", "_");
        }).style("background-color", function (d, i) {
            if (i === d.index) {
                return "rgba(0, 128, 0, 0.7)";
            } else if (d.column === "Potential threats" && d.detected === false && d.read_count < 10000) {
                return "rgba(255, 255, 255, 0)";
            } else if (d.column === "Detected" && d.detected === true) {
                return "rgba(255, 0, 0, 0.7)";
            }
        }).html(
            function (d) {
                if (d.column === "Not detected" && d.value === true) {
                    return d.species;
                } else if (d.column === "Potential threats" && d.detected === false && d.read_count < 10000) {
                    return d.species;
                } else if (d.column === "Detected" && d.detected === true) {
                    return d.species;
                }
            }
        );
    });
}