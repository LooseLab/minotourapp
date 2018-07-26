function getTotalReadsTable(flowCellId) {
    $.get("/table", {flowcellId: flowCellId}, result => {
        console.log(result);
        let data = JSON.parse(result.json);
        console.log(data);
        let table = $(".tableLand").DataTable({
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
    })
}