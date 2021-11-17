
class CnvChart {
  constructor (divId) {
    this._cnvChart = this._createCnvChart()
    this._divId = divId
  }

  get cnvChart () {
    return this._cnvChart
  }

  _createCnvChart () {
    return Highcharts.chart({
      chart: {
        type: `area`,
        renderTo: this._divId,
        reflow: true,
        marginLeft: 50,
        marginRight: 20,
        zoomType: `x`
      },
      boost: {
        useGPUTranslations: true,
        usePreAllocated: true
      }
    })
  }
}
