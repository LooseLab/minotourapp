class CnvChartController {
  constructor (divId, flowcellId) {
    this.cnvChart = new CnvChart(divId)
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

  _loadChartData (url, axiosInstance) {
    this.cnvChart.cnvChart.showLoading(`<div class="spinner-border text-success" role="status">
                        <span class = "sr-only"> Loading...</span></div>`)
    axiosInstance.get(url).then(
      response => {
        Object.entries(response.data).forEach(
          ([key, value]) => {
            if (key !== `plotting_data`) {
              this.cnvChart.cnvChart.addSeries({
                name: key,
                data: value,
                marker: {
                  symbol: `circle`,
                  radius: 2,
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
      }
    )
  }
}
