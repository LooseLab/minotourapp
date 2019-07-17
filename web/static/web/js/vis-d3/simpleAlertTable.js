function compare(a, b) {
// sorting function to put results in alphabetical order
    if (a["Validation species"] < b["Validation species"])
        return -1;
    if (a["Validation species"] > b["Validation species"])
        return 1;
    return 0;
}

function draw_simple_table(flowcellId) {
    // declare variables at scope start
    let table;
    let thead;
    let tbody;
    let rows;
    let cells;
    let columns = ["Status", "Low probability", "Validation species", "Detected"];
    let data;
    let limit;
    let flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');
    let barcode = get_selected_barcode();
    // if we are not on the metagenomics tab, don't waste time loading it
    if (flowcell_selected_tab_input.value !== "nav-metagenomics") {
        return;
    }
    // get the alerts 
    $.get("/api/v1/metagenomics/alerts", {flowcellId, barcode}, (result, statusText, xhr) => {
        // If there is no validation mapping results objects in the database for this run
        if (xhr.status == 204){
            d3.select("#validation").remove();
            return;
        }

        data = result.table;
        limit = result.conf_detect_limit;
        // sort the data aphabetically on the Validation species column
        data.sort(compare);
        // if the table has been created on the page drawing, select it so we can update the existing table
        if (d3.select(".simple-alert-table").classed("has-tabley?")) {
            table = d3.select(".simple-alert-table").select("table");
            thead = table.select("thead");
            tbody = table.select("tbody");
        }
        // else draw the table for the first time
        else {
            table = d3.select(".simple-alert-table").classed("has-tabley?", true)
                .style("width", "100%").append("table").attr("class", "table map-alert simple")
                .style("table-layout", "fixed").style("border-bottom", "1px solid #f4f4f4");
            thead = table.append('thead').append('tr');
            tbody = table.append('tbody').attr("class", "alert-tbody");
        }
        // select all the th elements
        thead.selectAll('th')
            // bind the columns array
            .data(columns).enter()
        // enter the selection, so enter is a selection of data that isn't already bound to an html element
            .append('th').style("text-align", "center")
            // give it an id matching the column name
            .attr("id", function (column) {
                return column.replace(" ", "_");
            }).attr("class", function (d, i) {
                // give it a class of it's index base 0 - used for to populate the table below
                return i;
                // make the status column 10% of the table width
            }).style("width", function(d){
                if (d==="Status"){
                    return "10%";
                }
                // return the column text
            }).text(function (column) {
                return column;
            });

        // select all the rows
        rows = tbody.selectAll('tr');
        // same binding of new elements as above
        rows
            .data(data)
            .enter()
            .append('tr').attr("id", function (d) {
                // give each row an id of the species title
            return d["Validation species"].replace(/ /g, "_");
        });
        // exit is a selection of html elements who's bound data in the array has been removed, so remove the html
        rows.data(data).exit().remove();
        // select all the new rows
        rows = tbody.selectAll('tr');
        // select all the cells in the rows
        cells = rows.selectAll('td');
        // for each column in each row, pass a dictionary down the chain of column and that columns value in each row
        cells.data(function (row) {
            return columns.map(function (column) {
                return {column: column, value: row[column]};
            });
        })
        // enter the new data, anything not bound to a html element append a cell
            .enter()
            .append('td');
        // select all the newly created and exisiting html td elements
        cells = rows.selectAll("td");
        // for each row return a mapped dictionary of all possible values for that rows columns
        cells.data(function (row) {
            return columns.map(function (column) {
                return {
                    column, value: row[column], read_count: row.read_count,
                    detected: row.Detected, species: row["Validation species"],
                    conf_limit: row.conf_limit, detected_at: row.detected_at
                };
            });
            // set the id of the cell as the column id
        }).attr("id", function (d, i) {
            return columns[i].replace(" ", "_");
        }).style("background-color", function (d) {
            // change status cell colour depending on findings
            if (d.column === "Status" && d.detected === false) {
                return "green";
            } else if (d.column === "Status" && d.detected === true) {
                // a nice pastely red colour
                return "rgba(255, 0, 0, 0.7)";
            }
        }).html(
            function (d) {
                // set the validation species html to not have a underscore
                d3.select("#Validation species").html("Validation species");
                // the code below places the html in the left cell if we have seen lots of reads but no matches
                if (d.column === "Low probability" && d.detected === false && d.conf_limit >= 50000) {
                    d3.select("#Low_probability").html("Low risk, less than 1 in " + limit + " reads");
                    return d.species;
                } else if (d.column === "Validation species" && d.detected === false && d.conf_limit < 50000){
                    d3.select("#Potential_threats").html("Validation species (< 1 in " + limit + " reads)");
                    return d.species;
                }
                else if (d.column === "Detected" && d.detected === true) {
                    return d.species + " (1 in "+ d.detected_at + " reads)";
                }
            }
        );
    });
}