/**
 * Controller for the tables on the Flowcell page and the Flowcell Manager page.
 */
class FlowcellTableController {
  /**
     *
     * @param linkDestination {string} The destination of the formatted link. flowcell_manager or flowcells.
     * @param elementId {string} The target Element to render this table inside.
     */
  constructor (linkDestination, elementId) {
    this._elementId = elementId
    this._linkDestination = linkDestination
    this._datatableObj = null
    this._interval = setInterval(this.refreshTable, 30000, elementId)
    // ON page unload clear the interval, so when user navigates away from table
    $(window).on(`unload`, function () {
      clearInterval(this._interval)
    })
  }

  /**
     * @function refresh the flowcell table.
     */
  refreshTable (_elementId) {
    $(`#${_elementId}`).DataTable().ajax.reload(null, false)
  }

  /**
     * @function Draw the flowcells index table.
     */
  renderTable () {
    const that = this
    var table_all_runs

    this._datatableObj = $(`#${this._elementId}`).DataTable({
      // callback for each created row, add styling class so cursor is pointer showing it's a link
      createdRow: function (row, data, index) {
        $(row).addClass(`stylable`)
        $(row).on(`click`, () => {
          // When we draw a row, make the whole row clickable, links to href destination
          window.location.href = `/web/private/${that._linkDestination}/${data.id}`
        })
      },
      // sScrollX: `100%`,
      // scrollX: true,
      order: [[3, `desc`]],
      ajax: {
        url: `/api/v1/flowcells/`,
        async: true,
        error: (xhr, error, code) => {
          console.error(xhr)
          console.error(code)
        }
      },
      columnDefs: [
        {
          targets: 0,
          data: `name`

        },
        {
          targets: 1,
          data: `size`
        },
        {
          targets: 2,
          data: `sample_name`
        },
        {
          targets: 3,
          data: `start_time`,
          render: function (data, type, full, meta) {
            return `${new Date(data).toISOString()}`
          }
        },
        {
          targets: 4,
          data: `number_reads`,
          render: $.fn.dataTable.render.number(`,`, `.`)
        },
        {
          targets: 5,
          data: `number_runs`
        },
        {
          targets: 6,
          data: `number_barcodes`
        },
        {
          targets: 7,
          data: `total_read_length`,
          render: function (data, type, row) {
            // The below function returns the yield in a human readable format
            if (type === `display` || type === `filter`) {
              var UNITS = [``, `k`, `M`, `G`, `T`, `P`, `E`, `Z`]
              var factor = 1000
              var suffix = `b`
              for (var i = 0; i < UNITS.length; i++) {
                if (Math.abs(data) < factor) {
                  data = data.toFixed(2)
                  if (/\.00$/.test(data)) {
                    data = data.substr(0, data.length - 3)
                  }
                  return `${`${data} ${UNITS[i]}${suffix}`}`
                }
                data /= factor
              }
              data = data.toFixed(2).replace(/\B(?=(\d{3})+(?!\d))/g, `,`)
              if (/\.00$/.test(data)) { data = data.substr(0, data.length - 3) }
              return `<a href="/web/private/${that._linkDestination}/${row.id}/">${`${data} Y${suffix}`}</a>`
            } else {
              return data
            }
          }
        },
        {
          targets: 8,
          data: `average_read_length`,
          render: $.fn.dataTable.render.number(`,`, `.`)

        },
        {
          targets: 9,
          data: `is_active`,
          render: (data, type, full, meta) => data.active ? `Active` : `Inactive`
        },
        {
          targets: 10,
          data: `archived`,
          render: (data, type, full, meta) => data.archived ? `Yes` : `No`

        },
        {
          targets: 11,
          data: `owner`,
          render: (data, type, full, meta) => data.owner ? `Yes` : `No`

        },
        {
          targets: 12,
          data: `permission`
        }
      ]
    })
  }
}
