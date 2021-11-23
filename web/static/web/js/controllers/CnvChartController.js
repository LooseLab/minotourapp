class CnvChartController {
  constructor (divId, flowcellId, detailDivId) {
    this.cnvChart = new CnvChart(divId)
    this._cnvDetailChart = new CnvChart(detailDivId)
    this._barcodeCnvSelect = $(`#barcodeCNVSelect`)
    this._barcodesLoaded = new Set()
    this._axiosInstance = axios.create({
      headers: { 'X-CSRFToken': getCookie(`csrftoken`) }
    })
    this._barcodeCnvSelect.on(
      `change`, (event) => {
        const barcodePk = this._barcodeCnvSelect.find(`:selected`).attr(`data-pk`)
        const url = `/api/v1/alignment/${flowcellId}/cnv-chart/${barcodePk}/2`
        this._cnvSelectChange(event, url, this._axiosInstance)
      })
    this._cnvDetailContigSelect = $(`#contigCNVDetail`)
    this._cnvDetailContigSelect.on(
      `change`, event => {
        const contigName = this._cnvDetailContigSelect.find(`:selected`).val()
        const barcodePk = this._barcodeCnvSelect.find(`:selected`).attr(`data-pk`)
        const minDiff = $(`#penalty-value`).val()
        const penValue = $(`#min-size-value`).val()
        const url = `/api/v1/alignment/${flowcellId}/cnv-chart-detail/${barcodePk}/${contigName}/${penValue}/${minDiff}`
        this._contigSelectChange(event, url, this._axiosInstance)
      }
    )
  }

  get cnvDetailChart () {
    return this._cnvDetailChart
  }

  _contigSelectChange (event, url, axiosInstance) {
    this._cnvDetailContigSelect.find(`#cnv-deet-placeholder`).remove()
    this._loadDetailChartData(url, axiosInstance)
  }

  _cnvSelectChange (event, url, axiosInstance) {
    this._barcodeCnvSelect.find(`#cnv-placeholder`).remove()
    this._loadChartData(url, axiosInstance)
  }

  _populateCnvBarcodeSelect (url, axiosInstance) {
    axiosInstance.get(url).then(
      response => {
        const barcodeHtml = []
        response.data.forEach(([barcode, barcodePk]) => {
          if (!this._barcodesLoaded.has(barcode)) {
            barcodeHtml.push(`<option data-pk="${barcodePk}">${barcode}</option>`)
          }
        })
        this._barcodeCnvSelect.append(barcodeHtml)
      })
  }

  _loadDetailChartData (url, axiosInstance) {
    this._cnvDetailChart.cnvChart.showLoading(`<div class="spinner-border text-success" role="status">
                        <span class = "sr-only"> Loading...</span></div>`)
    this._cnvDetailContigSelect.prop(`disabled`, true)
    if (this._cnvDetailChart.cnvChart.series.length) {
      this._cnvDetailChart.cnvChart.series[0].remove(false, false)
    }
    const colors = [`#f7c9c9`, `#B8e2f2`]

    axiosInstance.get(url).then(
      response => {
        console.log(response.data)
        const plotBands = response.data.plot_bands
        const points = response.data.points
        this._cnvDetailChart.cnvChart.addSeries({
          name: response.data.contig,
          data: points,
          color: `#000`,
          marker: {
            symbol: `circle`,
            radius: 2
          }
        }, false, false)
        this._cnvDetailChart.cnvChart.yAxis[0].addPlotLine({
          color: `black`,
          width: 2,
          value: 2,
          zIndex: 3
        })
        plotBands.forEach((band, ind) => {
          // console.log(band)
          // console.log(ind)
          // console.log(plotBands.length)
          // console.log(ind !== plotBands.length-1)
          if (ind !== plotBands.length - 1) {
            this._cnvDetailChart.cnvChart.xAxis[0].addPlotBand(
              {
                id: `band${ind}`,
                from: band,
                to: plotBands[ind + 1],
                color: colors[ind % 2],
                zIndex: 2
              }
            )
          }
        })
        this._cnvDetailChart.cnvChart.hideLoading()
        this._cnvDetailChart.cnvChart.redraw()
        this._cnvDetailContigSelect.prop(`disabled`, false)
      }
    )
  }

  _loadChartData (url, axiosInstance) {
    this.cnvChart.cnvChart.showLoading(`<div class="spinner-border text-success" role="status">
                        <span class = "sr-only"> Loading...</span></div>`)
    this._cnvDetailContigSelect.prop(`disabled`, true)
    while (this.cnvChart.cnvChart.length > 0) {
      this.cnvChart.cnvChart.series[0].remove(false, false)
    }

    axiosInstance.get(url).then(
      response => {
        const cnvContigDropDowns = [`<option id="cnv-deet-placeholder">Please choose...</option>`]
        Object.entries(response.data).forEach(
          ([key, value]) => {
            if (key !== `plotting_data`) {
              cnvContigDropDowns.push(`<option value="${key}">${key}</option>`)
              this.cnvChart.cnvChart.addSeries({
                name: key,
                data: value,
                marker: {
                  symbol: `circle`,
                  radius: 2
                }
              }, false, false)
            }
          })
        this.cnvChart.cnvChart.yAxis[0].addPlotLine({
          color: `black`,
          width: 2,
          value: 2,
          zIndex: 3
        })
        this.cnvChart.cnvChart.hideLoading()
        this.cnvChart.cnvChart.redraw()
        this._cnvDetailContigSelect.html(cnvContigDropDowns)
        this._cnvDetailContigSelect.prop(`disabled`, false)
      }
    )
  }
}
