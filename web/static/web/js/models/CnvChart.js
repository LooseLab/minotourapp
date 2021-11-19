
class CnvChart {
  constructor (divId) {
    this._divId = divId
    this._cnvChart = this._createCnvChart()
  }

  get cnvChart () {
    return this._cnvChart
  }

  _createCnvChart () {
    return Highcharts.chart({
      chart: {
        type: `scatter`,
        renderTo: this._divId,
        reflow: true,
        marginLeft: 50,
        marginRight: 20,
        zoomType: `xy`,
        events: {
          redraw () {
            this.series.forEach(singleSeries => {
              singleSeries.legendSymbol.attr({
                x: -2,
                y: 4,
                width: 15,
                height: 15
              })
            })
          }
        }
      },
      boost: {
        useGPUTranslations: true,
        usePreAllocated: true
      },
      plotOptions: {

      },
      title: {
        text: `CNV ME`
      },
      legend: {
        symbolHeight: 12,
        symbolWidth: 12,
        symbolRadius: 6
      }
    })
  }
}
