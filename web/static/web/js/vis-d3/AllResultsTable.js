function getTotalReadsTable(flowCellId) {
    $.get("/table", {flowcellId: flowCellId}, result => {
        console.log(result);
        let data = JSON.parse(result.json);
        let table = $(".tableLand");
        console.log(data);
        if (table.length !== 1) {
            table.DataTable().clear();
            table.DataTable().rows.add(data);
            table.DataTable().draw();
        } else {
            table.DataTable({
                    data: data, scrollY: "500px",
                    paging: false,
                    columns: [
                        {"data": "superkingdom"},
                        {"data": "phylum"},
                        {"data": "classy"},
                        {"data": "order"},
                        {"data": "family"},
                        {"data": "genus"},
                        {"data": "species"}
                    ]
                }
            )
        }

    })
}