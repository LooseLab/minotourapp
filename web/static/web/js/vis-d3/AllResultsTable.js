// function format(d) {
//
//     // `d` is the original data object for the row
//     return '<tr style="display: inline-block; width: 100%">' +
//         '<td style="width: 15px"></td>' +
//         '<td class="dropdown"> Kingdom (Num reads / Proportion <br>' + d.summed_superkingdom + ' / ' + d.prop_superkingdom + '%)    </td>' +
//         '<td class="dropdown"> Phylum (Num reads/ Proportion <br>' + d.summed_phylum + ' / ' + d.prop_phylum + '%)       </td>' +
//         '<td class="dropdown"> Class (Num reads / Proportion <br>' + d.summed_classy + ' / ' + d.prop_class + '%)   </td>' +
//         '<td class="dropdown"> Order (Num reads / Proportion <br>' + d.summed_order + ' / ' + d.prop_order + '%)    </td>' +
//         '<td class="dropdown"> Family (Num reads / Proprtion <br>' + d.summed_family + ' / ' + d.prop_family + '%)      </td>' +
//         '<td class="dropdown"> Genus (Num reads / Proportion <br>' + d.summed_genus + ' / ' + d.prop_genus + '%)     </td>' +
//         '<td class="dropdown"> Species (Num reads / Proportion <br>' + d.summed_species + ' / ' + d.prop_species + '%)        </td>' +
//         '</tr>'
//         ;
// }

function getTotalReadsTable(flowCellId) {
    // Get ttoal reads table updates the total reads table at te bottom of the page
    // Get data from the api
    $.get("/table", {flowcellId: flowCellId, visType: "table"}, result => {
        // Jquery selector
        let table = $(".tableLand");
        // If the results of the api query is undefined, there are no results in the DB, so retrun here
        if(result === undefined){
            return;
        }
        // If the table already exists, use the DataTable APi to update in place
        if (table.length !== 1) {
            table.DataTable().clear();
            table.DataTable().rows.add(result);
            table.DataTable().draw();
        } else {
            // else the databale must be initialised
            table = $(".tableLand").DataTable({
                    data: result, scrollY: "500px",
                    paging: false,
                    select: "single",
                    columns: [
                        {"data": "superkingdom"},
                        {"data": "phylum"},
                        {"data": "classy"},
                        {"data": "order"},
                        {"data": "family"},
                        {"data": "genus"},
                        {"data": "species"},
                        {"data": "num_matches"},
                        {"data": "prop_species"},
                    ]
                }
            );
        }
        // $('.tableLand tbody').on('click', 'td.details-control', function () {
        //     var tr = $(this).closest('tr');
        //     var tdi = tr.find("i.fa");
        //     var row = table.row(tr);
        //
        //     if (row.child.isShown()) {
        //         // This row is already open - close it
        //         row.child.hide();
        //         tr.removeClass('shown');
        //         tdi.first().removeClass('fa-minus-square');
        //         tdi.first().addClass('fa-plus-square');
        //     }
        //     else {
        //         // Open this row
        //         row.child(format(row.data())).show();
        //         tr.addClass('shown');
        //         tdi.first().removeClass('fa-plus-square');
        //         tdi.first().addClass('fa-minus-square');
        //     }
        // });
        //
        // table.on("user-select", function (e, dt, type, cell, originalEvent) {
        //     if ($(cell.node()).hasClass("details-control")) {
        //         e.preventDefault();
        //     }
        // });
    });

}