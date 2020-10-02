class SummaryTabController {
  constructor (flowcellId) {
    this._flowcellId = flowcellId
    this._dataTableObj = null
    this._axiosInstance = axios.create({
      headers: { 'X-CSRFToken': getCookie(`csrftoken`) }
    })
  }

  /**
   * Update tab function, called from the FlowcellTabController
   */
  updateTab () {
    const htmlTables = [[`flowcell-run-basecalled-summary-html`, `#flowcell-run-basecalled-summary-html-div`], [`run-summaries-html`, `#run-summaries-div`]]
    this._dataTableObj ? this._dataTableObj.ajax.reload(null, false) : this._drawMinknowMessagesTable(this._flowcellId)
    htmlTables.forEach(([urlEnd, divId]) => {
      this._requestHtmlSummaries(this._flowcellId, urlEnd, divId)
    })
  }

  /**
   * Draw the minKnow messages data table on the flowcell summary tab
   * @private
   * @param flowcellId {number} The primary key of this flowcell
   */
  _drawMinknowMessagesTable (flowcellId) {
    this._dataTableObj = $(`#minknow-messages-table`).DataTable({
      initComplete: (settings, json) => {
        $(`#minknow-messages-tdiv`).css({ display: `block` })
        $(`#minknow-messages-div`).remove()
      },
      ajax: {
        url: `/api/v1/flowcells/${flowcellId}/minknow-messages/`,
        async: true,
        error: (xhr, error, code) => {
          console.error(xhr)
          console.error(code)
        }
      },
      order: [[1, `desc`]],
      columnDefs: [
        {
          targets: 0,
          data: `messageText`

        },
        {
          targets: 1,
          data: `timestamp`
        }]
    })
  }

  /**
   * Request html tables to be displayed on the summary tab, append them to specified Div
   * @param flowcellId {number} The primary key of this flowcell
   * @param urlEnd {string} The end of the url to request
   * @param divId {string} The id of the Div that we will append the html table to
   * @private
   */
  _requestHtmlSummaries (flowcellId, urlEnd, divId) {
    const runSummariesTable = $(`${divId}`)
    this._axiosInstance.get(`/api/v1/reads/flowcells/${flowcellId}/${urlEnd}/`, {
      responseType: `html`
    })
      .then(response => {
        runSummariesTable.removeClass()
        runSummariesTable.html(response.data)
      }
      ).catch(error => {
        console.error(`Error: ${error}`)
      })
  }
}
