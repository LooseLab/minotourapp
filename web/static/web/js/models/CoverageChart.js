class CoverageChart {
  /*
      Create a new coverage chart, which is a highcharts stepped line chart
      returns a getter, and defines the master chart detail chart after selection event
       */

  constructor (divId) {
    // the chart
    this._divIdMaster = `${divId}_master`
    this._divIdDetail = `${divId}_detail`

    this._masterChart = this.createMasterChart()
    this._detailChart = this.createDetailChart()
  }

  get masterChart () {
    return this._masterChart
  }

  get detailChart () {
    return this._detailChart
  }

  createMasterChart () {
    const self = this

    return Highcharts.chart({

      chart: {
        renderTo: this._divIdMaster,
        reflow: true,
        marginLeft: 50,
        marginRight: 20,
        zoomType: `x`,
        step: true,
        events: {
          selection: self.afterSelection.bind(self)
        }
      },
      title: {
        text: null
      },
      xAxis: {
        type: `line`,
        showLastTickLabel: true,
        plotBands: [{
          id: `mask-before`,
          color: `rgba(0, 0, 0, 0.2)`
        }],
        title: {
          text: null
        }
      },
      yAxis: {
        title: {
          text: null
        },
        visible: false
      },
      legend: {
        enabled: false
      },
      credits: {
        enabled: false
      },
      plotOptions: {
        series: {
          fillColor: {
            // linearGradient: [0, 0, 0, 70],
            stops: [
              [0, Highcharts.getOptions().colors[0]],
              [1, `rgba(255,255,255,0)`]
            ]
          },
          lineWidth: 1,
          marker: {
            enabled: false
          },
          shadow: false,
          states: {
            hover: {
              lineWidth: 1
            }
          },
          enableMouseTracking: false
        }
      },
      series: [{
        type: `area`,
        name: null,
        data: [],
        step: true
      }],
      exporting: {
        enabled: false
      }
    })
  }

  createDetailChart () {
    const self = this
    return Highcharts.chart({
      chart: {
        renderTo: this._divIdDetail,
        reflow: true,
        marginLeft: 50,
        marginRight: 20,
        style: {
          position: `absolute`
        },
        boost: {
          useGPUTranslations: true,
          usePreAllocated: true
        },
        resetZoomButton: {
          fill: `white`,
          stroke: `silver`,
          position: {
            // align: 'right', // by default
            // verticalAlign: 'top', // by default
            x: -10,
            y: 10
          },
          relativeTo: `chart`,
          r: 0,
          states: {
            hover: {
              fill: `#41739D`,
              style: {
                color: `white`
              }
            }
          }
        },
        step: true,
        zoomType: `x`,
        events: {
          selection: self.afterSelection.bind(self)
        }
      },
      credits: {
        enabled: false
      },
      title: {
        text: null
      },
      xAxis: {
        minRange: 1500,
        crosshair: true
      },
      yAxis: {
        crosshair: true,
        title: {
          text: null
        },
        min: 0
      },
      legend: {
        enabled: false
      },
      plotOptions: {
        series: {
          marker: {
            enabled: false,
            states: {
              hover: {
                enabled: true,
                radius: 3
              }
            }
          }
        }
      },
      series: [{
        type: `area`,
        name: null,
        data: [],
        step: `left`
      }],
      exporting: {
        enabled: true
      },
      tooltip: {
        pointFormat: `<span style="color:{point.color}">●</span> The coverage for bin starting {point.x}: <b>{point.y}</b><br/>`
      }
    })
  }

  /**
     * After selection event callback called when a zoom band is chosen on either the detail or master chart.
     * @param event {Object} Dom event object
     */
  afterSelection (event) {
    const min = Math.trunc(event.xAxis[0].min)
    const max = Math.trunc(event.xAxis[0].max)
    event.preventDefault()
    this._masterChart.xAxis[0].removePlotBand(`mask-before`)
    this._masterChart.xAxis[0].addPlotBand({
      id: `mask-before`,
      from: min,
      to: max,
      color: `rgba(0, 0, 0, 0.2)`
    })
    this._detailChart.xAxis[0].setExtremes(min, max)
  }
}
