/* global getCookie, CoverageChart, axios */
class CoverageChartController {
  // controller for coverage chart
  constructor (divId) {
    // constructs new Coverage chart class, and link to chromosome select field
    this._coverageChart = new CoverageChart(divId)
    this._oldSumToCheck = 0
    this._axiosInstance = axios.create({
      headers: { 'X-CSRFToken': getCookie(`csrftoken`) }
    })
  }

  /**
   * Reload the master and detail chart data.
   * @param url {string} The url to load data from
   */
  reloadCoverageCharts (url) {
    // reload the data on click to reset or chromosome change
    this.loadChartData(url)
  }

  /**
   * Make the ajax call to get chart data and insert it if it is not identical
   * @param url
   */
  loadChartData (url) {
    const self = this
    self._coverageChart.masterChart.showLoading(`<div class="spinner-border text-success" role="status">
                        <span class = "sr-only"> Loading...</span></div>`)

    this._axiosInstance.get(url).then(
      response => {
        const data = response.data.newChartData
        const sumToCheck = response.data.sumToCheck
        self._refLength = response.data.refLength
        if (this._oldSumToCheck !== sumToCheck) {
          console.log(`new data`)
          self._coverageChart.masterChart.xAxis[0].setExtremes(0, response.data.refLength)
          self._coverageChart.masterChart.series[0].setData(data.sequenced, false, false, false)
          self._coverageChart.masterChart.series[1].setData(data.unblocked, false, false, false)
          // self._coverageChart.detailChart.series[0].setData(data.Sequenced, false, false, false)
          // self._coverageChart.detailChart.series[1].setData(data.Unblocked, false, false, false)
          self._coverageChart.masterChart.xAxis[0].removePlotBand(`mask-before`)
          // self._coverageChart.detailChart.hideLoading()
          self._coverageChart.masterChart.hideLoading()
          // self._coverageChart.detailChart.redraw()
          self._coverageChart.masterChart.redraw()
        } else {
          // self._coverageChart.detailChart.hideLoading()
          self._coverageChart.masterChart.hideLoading()
        }
        this._oldSumToCheck = sumToCheck
      }).catch(
      error => {
        console.error(error)
      })
  }

  /**
   * Reset the chart zoom on the detail chart.
   */
  resetDetailChartZoom () {
    this._coverageChart.detailChart.xAxis[0].setExtremes(0, self._refLength)
    this._coverageChart.detailChart.series[0].setData(0,false, false, false)
    this._coverageChart.masterChart.xAxis[0].removePlotBand(`mask-before`)
  }
}
