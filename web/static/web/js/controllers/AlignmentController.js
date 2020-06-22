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
    this._selects = [this._readTypeSelect, this._referenceSelect, this._chromosomeSelect, this._barcodeSelect]
    this._addChangeListenerToSelects()
    this.requestMappedChromosomes(flowcellId)
    this._drawPafSummaryTable(flowcellId)
    this._requestColumnChartData(flowcellId)
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
      that._requestColumnChartData(flowcellId)
    }
  }

  /**
   *
   * @param flowcellId
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
      this.selectOnChange({ srcElement: { id: `referenceSelect` } })
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
        this.selectOnChange(event)
      })
    })
  }

  /**
   * Request the mapped data and add it to the the first
   * @param {number} flowcellId The primary key of the flowcell record in the database
   */
  requestMappedChromosomes (flowcellId) {
    this._axiosInstance.get(`/api/v1/alignment/${flowcellId}/mapped-references`).then(response => {
      const readTypes = new Set(response.data.map(el => {
        return `<option value="${el.referenceId}" data-jm-id="${el.jmId}">${el.referenceName}</option>`
      }))
      this._selectData = response.data
      console.log(readTypes)
      this._referenceSelect.append([...readTypes])
      // todo add code to append new options to the dropdowns here
    })
  }

  /**
   * Request the chromosomes that have reads mapped to using minimap2
   * and update the select box on tab Mapping
   */
  selectOnChange (event) {
    let url
    let barcodeName
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
      barcodeName = this._barcodeSelect.val()
      readTypeId = this._readTypeSelect.val()
      chromosomeId = this._chromosomeSelect.val()
      url = `/api/v1/alignment/coverage/${taskId}/${barcodeName}/${readTypeId}/${chromosomeId}`
      this.coverageChartController.reloadCoverageCharts(url)
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
            url: `/api/v1/alignment/${flowcellId}/pafsummarytable`,
            async: true,
            error: (xhr, error, code) => {
              console.error(xhr)
              console.error(code)
            }
          },
          columns: [
            { data: `barcode_name` },
            { data: `chromosome__reference__name` },
            { data: `chromosome__line_name` },
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

  _requestColumnChartData (flowcellId) {
    // TODO BELOW IS A HUGE NONO
    if (!this.chart_per_chrom_cov) {
      this.chart_per_chrom_cov = makeColumnChart(
        `per-chrom-cov`,
        `CHROMOSOME COVERAGE`,
        `CHROMOSOME COVERAGE`
      )
    }
    if (!this.chart_per_chrom_avg) {
      this.chart_per_chrom_avg = makeColumnChart(
        `per-chrom-avg`,
        `MEAN READ LENGTH BY CHROMOSOME`,
        `MEAN READ LENGTH BY CHROMOSOME`
      )
    }
    this._axiosInstance.get(`/api/v1/alignment/${flowcellId}/pafsummary/`).then(
      response => {
        const data = response.data
        const charts = [[this.chart_per_chrom_cov, `coverageData`], [this.chart_per_chrom_avg, `avgRLData`]]
        charts.forEach(([chart, key]) => {
          this.updateCoverageColumnCharts(chart, data[key], data.categories)
        })
      }
    ).catch(
      error => {
        console.error(error)
      })
  }

  /**
   * Update the coverage column chart data
   * @param chart {Object} The highcharts class representing the column chart we are updating
   * @param seriesData {Array.Object} Array of new series to insert into this chart
   * @param categories {Array.String} Array of category names
   */
  updateCoverageColumnCharts (chart, seriesData, categories) {
    console.log(chart)
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
