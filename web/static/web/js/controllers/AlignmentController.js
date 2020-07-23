/* global CoverageChartController, getCookie, axios */
class AlignmentController {
  /**
   * Controller for the mapping tab.
   * @param flowcellId {number} Primary key of flowcell record in database.
   */
  constructor (flowcellId) {
    this._flowcellId = flowcellId
    this.coverageChartController = new CoverageChartController(`coverage_div`, ``)
    this._axiosInstance = axios.create({
      headers: { 'X-CSRFToken': getCookie(`csrftoken`) }
    })
    this._selectData = []
    this._readTypeSelect = $(`#readTypeSelect`)
    this._referenceSelect = $(`#referenceSelect`)
    this._chromosomeSelect = $(`#chromosomeSelect`)
    this._barcodeSelect = $(`#barcodeSelect`)
    this._lookupDownStreamSelects = {
      referenceSelect: [this._readTypeSelect, this._chromosomeSelect, this._barcodeSelect],
      readTypeSelect: [this._chromosomeSelect, this._barcodeSelect],
      chromosomeSelect: [this._barcodeSelect],
      barcodeSelect: []
    }
    this._sumToCheck = 0
    this._selects = [this._readTypeSelect, this._referenceSelect, this._chromosomeSelect, this._barcodeSelect]
    this._addChangeListenerToSelects()
    // this._requestMappedChromosomes(flowcellId)
    // this._drawPafSummaryTable(flowcellId)
    this._createColumnCharts(flowcellId)
    this._addListenerToResetButton()
    this._interval = setInterval(this._reloadPageData, 45000, flowcellId, this)
    $(window).on(`unload`, function () {
      console.log(`clearing Alignment interval`)
      clearInterval(this._interval)
    })
  }

  /**
   * Called from flowcell tab controller when we switch to that tab
   */
  updateTab () {
    this._reloadPageData(this._flowcellId, this)
  }

  /**
   * Reload the page data on interval
   * @param flowcellId {number} The primary key of the flowcellID
   * @param that {AlignmentController} Reference to alignment controllers scope.
   * @private
   */
  _reloadPageData (flowcellId, that) {
    if (getSelectedTab() === `sequence-mapping`) {
      that._updateOptionDropDowns(flowcellId)
      that._drawPafSummaryTable(flowcellId)
      that._createColumnCharts(flowcellId)
      that._fetchChromosomeInGenomeCoverageData(flowcellId)
    }
  }

  /**
   * Update the drop downs, setting values for uniqueness then adding them to class
   * @param flowcellId {number} The flowcell primary key record
   * @private
   */
  _updateOptionDropDowns (flowcellId) {
    this._axiosInstance.get(`/api/v1/alignment/${flowcellId}/mapped-references`).then(response => {
      const referenceId = $(`#referenceSelect`).val()
      const readTypes = new Set(response.data.map(el => {
        if (el.referenceId === parseInt(referenceId)) {
          return `<option value="${el.referenceId}" data-jm-id="${el.jmId}" selected>${el.referenceName}</option>`
        } else {
          return `<option value="${el.referenceId}" data-jm-id="${el.jmId}">${el.referenceName}</option>`
        }
      }))
      this._selectData = response.data
      this._referenceSelect.html([...readTypes])
      this._selectOnChange({ srcElement: { id: `referenceSelect` } }, false)
    }).catch(error => {
      console.error(error)
    })
  }

  /**
   * Add the zoom reset to the rest button
   * @private
   */
  _addListenerToResetButton () {
    $(`#reset-chart-zoom`).on(`click`, () => {
      this.coverageChartController.resetDetailChartZoom()
    })
  }

  /**
   * Add on change listener to the select dropDowns
   * @private
   */
  _addChangeListenerToSelects () {
    this._selects.forEach(select => {
      select.on(`change`, () => {
        this._selectOnChange(event, true)
      })
    })
  }

  /**
   * Request the chromosomes that have reads mapped to using minimap2
   * and update the select box on tab Mapping
   * @param userActivated {boolean} If the user has changed the select, not the interval call
   */
  _selectOnChange (event, userActivated) {
    let url
    let barcodeId
    let readTypeId
    let chromosomeId
    const changedSelectId = event.srcElement.id
    const taskId = $(`#referenceSelect option:selected`).attr(`data-jm-id`)
    // Filter the data to only be for the selected value downstream via jobmaster pk
    const filteredData = this._selectData.filter(el => {
      return el.jmId === parseInt(taskId)
    })
    const downStreamers = this._lookupDownStreamSelects[changedSelectId]
    // remove the choose placeholder option
    $(`#${changedSelectId.substring(0, changedSelectId.length - 6)}Placeholder`).remove()
    downStreamers.forEach(select => {
      const currentSelection = select.val()
      const selectId = select.attr(`id`)
      const base = selectId.substring(0, selectId.length - 6)
      const id = `${base}Id`
      const name = `${base}Name`
      let appropriateOptionExists = false
      const options = new Set(filteredData.map(el => {
        if (el[id] === parseInt(currentSelection) || (base === `barcode` && el[id] === currentSelection)) {
          appropriateOptionExists = true
          return `<option value="${el[id]}" selected>${el[name]}</option>`
        } else {
          return `<option value="${el[id]}">${el[name]}</option>`
        }
      }))
      if (!appropriateOptionExists) {
        options.add(`<option value="-1" id="${base}Placeholder" selected>Please Choose</option>`)
      }
      select.html([...options].reverse())
      // if no selection has been made, we have no task id
    })
    // if we have children in every select and the selected option isn't a please choose
    if (this._selects.every(select => {
      return Boolean(select.children().length) === true && select.val() !== `-1`
    })) {
      // undisable the rest button
      $(`#reset-chart-zoom`).attr(`disabled`, false)
      barcodeId = this._barcodeSelect.val()
      readTypeId = this._readTypeSelect.val()
      chromosomeId = this._chromosomeSelect.val()
      url = `/api/v1/alignment/coverage/${taskId}/${barcodeId}/${readTypeId}/${chromosomeId}`
      if (userActivated) {
        this.coverageChartController.reloadCoverageCharts(url)
        this.coverageChartController.resetDetailChartZoom()
      }
      this._fetchBarcodeCoverageColumnChartsData(this._flowcellId, chromosomeId, true, barcodeId)
    }
  }

  /**
   * Draw or update the PafSummary table on the mapping tab.
   * @param flowcellId {number} Pk of the flowcell record in the database.
   * @private
   */
  _drawPafSummaryTable (flowcellId) {
    // Get total reads table updates the total reads table at te bottom of the page
    // Get data from the api
    // Jquery selector
    const table = $(`.mapping-table`)
    // If the table already exists, use the DataTable APi to update in place
    if ($.fn.DataTable.isDataTable(table)) {
      table.DataTable().ajax.reload(null, false)
    } else {
      // else the databale must be initialised
      table.DataTable(
        {
          ajax: {
            url: `/api/v1/alignment/${flowcellId}/paf-summary-table`,
            async: true,
            error: (xhr, error, code) => {
              console.error(xhr)
              console.error(code)
            }
          },
          columns: [
            { data: `barcode_name` },
            { data: `reference_name` },
            { data: `chromosome_name` },
            { data: `reference_line_length` },
            { data: `read_count` },
            { data: `total_yield` },
            { data: `average_read_length` },
            { data: `coverage` }
          ]
        }
      )
    }
  }

  _createColumnCharts () {
    const options = { lang: { noData: `Please select above` } }
    if (!this.chart_per_barcode_cov) {
      this.chart_per_barcode_cov = makeColumnChart(
        `per-barcode-cov`,
        `BARCODE COVERAGE`,
        `BARCODE COVERAGE`
      )
    }
    if (!this.chart_per_barcode_avg) {
      this.chart_per_barcode_avg = makeColumnChart(
        `per-barcode-avg`,
        `MEAN READ LENGTH BY BARCODE`,
        `MEAN READ LENGTH BY BARCODE`
      )
    }
    if (!this.chart_per_genome_cov) {
      this.chart_per_genome_cov = makeColumnChart(
        `per-genome-cov`,
        `GENOME COVERAGE`,
        `GENOME COVERAGE`
      )
    }
    if (!this.chart_per_genome_avg) {
      this.chart_per_genome_avg = makeColumnChart(
        `per-genome-avg`,
        `MEAN READ LENGTH BY GENOME`,
        `MEAN READ LENGTH BY GENOME`
      )
    }
    this.chart_per_barcode_avg.update(options)
    this.chart_per_barcode_cov.update(options)
  }

  /**
   * Fetch the data for the per barcode in a chromosome chart data
   * @param flowcellId {number} The primary key of the flowcell Record
   * @param chromosomeId {number} The primary key of the Chromosome Record
   * @param changedBarcode {boolean} Redraw is result of user changing barcode
   * @param barcodePk {number} The primary key of the barcode
   * @private
   */
  _fetchBarcodeCoverageColumnChartsData (flowcellId, chromosomeId, changedBarcode, barcodePk) {
    this._axiosInstance.get(`/api/v1/alignment/${flowcellId}/pafsummary/`, { params: { chromosomeId, barcodePk } }).then(
      response => {
        const data = response.data
        const charts = [[this.chart_per_barcode_cov, `coverageData`], [this.chart_per_barcode_avg, `avgRLData`]]
        charts.forEach(([chart, key]) => {
          this._updateCoverageColumnCharts(chart, data[key], data.categories, changedBarcode)
        })
      }
    ).catch(
      error => {
        console.error(error)
      })
  }

  /** TODO almost duplicate of above, merge?
   * TODO merge charts into one grouped chart
   * Fetch per chromosome in Genome data for a flowcell with mapping tasks
   * @param flowcellId {number} The primary key of the flowcell record
   * @private
   */
  _fetchChromosomeInGenomeCoverageData (flowcellId) {
    this._axiosInstance.get(`/api/v1/alignment/${flowcellId}/genome-coverage-summary/`).then(
      response => {
        const data = response.data
        const charts = [[this.chart_per_genome_cov, `coverageData`], [this.chart_per_genome_avg, `avgRLData`]]
        charts.forEach(([chart, key]) => {
          this._updateCoverageColumnCharts(chart, data[key], data.categories, false)
        })
      }
    )
  }

  /**
   * Update the column chart data
   * @param chart {Object} The highcharts class representing the column chart we are updating
   * @param seriesData {Array.Object} Array of new series to insert into this chart
   * @param categories {Array.String} Array of category names
   * @param changedBarcode {boolean} If redraw is a result of barcode change
   */
  _updateCoverageColumnCharts (chart, seriesData, categories, changedBarcode) {
    if (changedBarcode){
      while (chart.series.length) {
        chart.series[0].remove(false)
      }
    }
    chart.xAxis[0].setCategories(categories)
    seriesData.forEach(series => {
      const oldSeries = chart.series.filter(el => {
        return el.name === series.name
      })
      if (oldSeries.length) {
        if (!checkHighChartsDataIsIdentical(series.data, oldSeries[0].options.data)) {
          oldSeries[0].setData(series.data)
        }
      } else {
        chart.addSeries(series)
      }
    })
  }
}
