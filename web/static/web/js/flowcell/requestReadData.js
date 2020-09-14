function requestReadData (id) {
  var table

  if ($.fn.dataTable.isDataTable(`#example`)) {
    table = $(`#example`).DataTable()
  } else {
    table = $(`#example`).DataTable({
      processing: true,
      serverSide: true,
      pagingType: `simple`,
      deferLoading: 0,
      ajax: {
        url: `/web/private/flowcells_reads_data/`,
        data: {
          flowcell_id: id
        }
      },
      columns: [
        { data: `read_id` },
        { data: `read` },
        { data: `channel` },
        { data: `sequence_length` },
        { data: `run__runid` },
        { data: `barcode__name` },
        {
          className: `details-control`,
          sortable: false,
          data: null,
          defaultContent: `<button class="btn btn-primary" type="button">Show sequence</button>`
        }
      ]
    })

    $(`#example`).on(`click`, `td.details-control`, function () {
      var tr = $(this).closest(`tr`)
      var row = table.row(tr)

      if (row.child.isShown()) {
        row.child.hide()
        tr.removeClass(`shown`)
      } else {
        row.child(format(row.data())).show()
        tr.addClass(`shown`)
      }
    })
  }
}
