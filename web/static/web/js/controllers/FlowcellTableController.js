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
    console.log(`[FlowcellTable] constructor`, { linkDestination, elementId, ts: new Date().toISOString() })
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
    console.log(`[FlowcellTable] refresh tick`, { elementId: _elementId, ts: new Date().toISOString() })
    $(`#${_elementId}`).DataTable().ajax.reload(null, false)
  }

  /**
     * @function Draw the flowcells index table.
     */
  renderTable () {
    const that = this
    const messageEl = $(`#messages`)
    console.log(`[FlowcellTable] renderTable start`, { elementId: this._elementId, ts: new Date().toISOString() })
    this._datatableObj = $(`#${this._elementId}`).DataTable({
      processing: true,
      // callback for each created row, add styling class so cursor is pointer showing it's a link
      createdRow: (row, data, index) => {
        $(row).addClass(`stylable`)
        $(row).on(`click`, () => {
          // When we draw a row, make the whole row clickable, links to href destination
          window.location.href = `/web/private/${that._linkDestination}/${data.id}`
        })
      },
      responsive: {
        details: false
      },
      order: [[3, `desc`]],
      ajax: {
        url: `/api/v1/reads/flowcells/`,
        async: true,
        dataSrc: (json) => {
          console.log(`[FlowcellTable] ajax dataSrc`, { hasJson: !!json, keys: json ? Object.keys(json) : [], ts: new Date().toISOString() })
          if (!json || !Array.isArray(json.data)) {
            console.error(`Unexpected flowcells payload`, json)
            if (messageEl.length) {
              messageEl.html(`<div class="alert alert-danger">Flowcell API returned an unexpected payload format.</div>`)
            }
            return []
          }
          console.log(`[FlowcellTable] ajax rows`, { count: json.data.length, firstRow: json.data[0], ts: new Date().toISOString() })
          return json.data
        },
        error: (xhr, error, code) => {
          console.error(xhr)
          console.error(code)
          console.error(`[FlowcellTable] ajax error`, { status: xhr ? xhr.status : null, error, code, responseText: xhr ? xhr.responseText : null, ts: new Date().toISOString() })
          if (messageEl.length) {
            const status = xhr && xhr.status ? xhr.status : `unknown`
            messageEl.html(`<div class="alert alert-danger">Failed to load flowcells (HTTP ${status}). Check browser console/network for details.</div>`)
          }
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
          render: (data, type, full, meta) => {
            const dt = new Date(data)
            if (Number.isNaN(dt.getTime())) {
              return data || ``
            }
            return dt.toISOString()
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
          render: (data, type, row) => {
            return humanReadableYield(data, type, row)
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
          render: (data, type, full, meta) => {
            return data ? `Active` : `Inactive`
          }
        },
        {
          targets: 10,
          data: `archived`,
          render: (data, type, full, meta) => {
            return data ? `Yes` : `No`
          }

        },
        {
          targets: 11,
          data: `owner`,
          render: (data, type, full, meta) => {
            return data ? `Yes` : `No`
          }

        },
        {
          targets: 12,
          data: `permission`
        }
      ]
    })
    this._datatableObj.on(`init.dt`, () => {
      console.log(`[FlowcellTable] init.dt`, { rows: this._datatableObj.rows().count(), ts: new Date().toISOString() })
    })
    this._datatableObj.on(`xhr.dt`, (e, settings, json, xhr) => {
      console.log(`[FlowcellTable] xhr.dt`, {
        status: xhr ? xhr.status : null,
        rows: json && Array.isArray(json.data) ? json.data.length : null,
        ts: new Date().toISOString()
      })
    })
    this._datatableObj.on(`draw.dt`, () => {
      console.log(`[FlowcellTable] draw.dt`, { rows: this._datatableObj.rows().count(), ts: new Date().toISOString() })
    })
  }
}
