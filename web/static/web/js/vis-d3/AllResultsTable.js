function format(d) {

    // `d` is the original data object for the row
    return '<tr style="display: inline-block; width: 100%">' +
        '<td style="width: 15px"></td>' +
        '<td class="dropdown"> Kingdom Num reads / Proportion <br>' + d.summed_superkingdom + ' / ' + d.prop_superkingdom + '%    </td>' +
        '<td class="dropdown"> Phylum Num reads/ Proportion <br>' + d.summed_phylum + ' / ' + d.prop_phylum + '%       </td>' +
        '<td class="dropdown"> Class Num reads / Proportion <br>' + d.summed_classy + ' / ' + d.prop_class + '%   </td>' +
        '<td class="dropdown"> Order Num reads / Proportion <br>' + d.summed_order + ' / ' + d.prop_order + '%    </td>' +
        '<td class="dropdown"> Family Num reads / Proprtion <br>' + d.summed_family + ' / ' + d.prop_family + '%      </td>' +
        '<td class="dropdown"> Genus Num reads / Proportion <br>' + d.summed_genus + ' / ' + d.prop_genus + '%     </td>' +
        '<td class="dropdown"> Species Num reads / Proportion <br>' + d.summed_species + ' / ' + d.prop_species + '%        </td>' +
        '</tr>'
;}

function getTotalReadsTable(flowCellId) {
    $.get("/table", {flowcellId: flowCellId, visType: "table"}, result => {
        console.log(result);
        let data = JSON.parse(result.json);
        let table = $(".tableLand");

        console.log(data);
        if (table.length !== 1) {
            table.DataTable().clear();
            table.DataTable().rows.add(data);
            table.DataTable().draw();
        } else {
            table = $(".tableLand").DataTable({
                    data: data, scrollY: "500px",
                    paging: false,
                select:"single",
                    columns: [
                        {
                            "className": 'details-control',
                            "orderable": false,
                            "data": null,
                            "defaultContent": '',
                            "render": function () {
                                return '<i class="fa fa-plus-square" aria-hidden="true"></i>';
                            },
                            width: "15px"
                        },
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
    $('.tableLand tbody').on('click', 'td.details-control', function () {
        console.log("clicked");
        var tr = $(this).closest('tr');
             var tdi = tr.find("i.fa");
             var row = table.row(tr);

             if (row.child.isShown()) {
                 // This row is already open - close it
                 console.log("closing row");
                 row.child.hide();
                 tr.removeClass('shown');
                 tdi.first().removeClass('fa-minus-square');
                 tdi.first().addClass('fa-plus-square');
             }
             else {
                 // Open this row
                 console.log("opening row");
                 row.child(format(row.data())).show();
                 tr.addClass('shown');
                 tdi.first().removeClass('fa-plus-square');
                 tdi.first().addClass('fa-minus-square');
             }
         });

         table.on("user-select", function (e, dt, type, cell, originalEvent) {
             if ($(cell.node()).hasClass("details-control")) {
                 e.preventDefault();
             }
         });
    })

}