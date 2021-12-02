class CoverageChart {
  /*
      Create a new coverage chart, which is a highcharts stepped line chart
      returns a getter, and defines the master chart detail chart after selection event
       */

  constructor (divId) {
    // the chart
    this._divIdMaster = `${divId}_master`
    this._divIdDetail = `${divId}_detail`
    this._axiosInstance = axios.create({
      headers: { 'X-CSRFToken': getCookie(`csrftoken`) }
    })
    this._masterChart = this.createMasterChart()
    this._detailChart = this.createDetailChart()
    this._readTypeSelect = $(`#readTypeSelect`)
    this._referenceSelect = $(`#referenceSelect`)
    this._chromosomeSelect = $(`#chromosomeSelect`)
    this._barcodeSelect = $(`#barcodeSelect`)
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
        type: `area`,
        renderTo: this._divIdMaster,
        reflow: true,
        marginLeft: 50,
        marginRight: 20,
        zoomType: `x`,
        events: {
          selection: self.afterSelection.bind(self)
        }
      },
      title: {
        text: null
      },
      xAxis: {
        type: `area`,
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
      boost: {
        useGPUTranslations: true,
        usePreAllocated: true
      },
      plotOptions: {
        series: {
          stacking: `normal`,
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
        name: `sequenced`,
        data: [],
        step: true
      }, {
        name: `unblocked`,
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
        type: `area`,
        renderTo: this._divIdDetail,
        reflow: true,
        marginLeft: 50,
        marginRight: 20,

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
        minRange: 100,
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
        enabled: true
      },
      plotOptions: {
        series: {
          stacking: `normal`,
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
        name: `sequenced`,
        data: [],
        step: true
      }, {
        name: `unblocked`,
        data: [],
        step: true
      }],
      exporting: {
        enabled: false
      },
      tooltip: {
        pointFormatter: function () {
          return `<span style="color:{this.color}">●</span> The coverage for bin ${this.x} - ${this.x + 10}: <b>${this.y}</b><br/>`
        }
      }
    })
  }

  /**
     * After selection event callback called when a zoom band is chosen on either the detail or master chart.
     * @param event {Object} Dom event object
     */
  afterSelection (event) {
    event.preventDefault()
    const self = this
    console.log(event.xAxis)
    let min = Math.trunc(event.xAxis[0].min)
    const max = Math.trunc(event.xAxis[0].max)
    const taskId = $(`#referenceSelect option:selected`).attr(`data-jm-id`)
    const barcodeId = this._barcodeSelect.val()
    const readTypeId = this._readTypeSelect.val()
    const chromosomeId = this._chromosomeSelect.val()
    const url = `/api/v1/alignment/coverage/${taskId}/${barcodeId}/${readTypeId}/${chromosomeId}/${min}/${max}`
    min = min < 0 ? 0 : min
    console.log(min)
    this._masterChart.xAxis[0].removePlotBand(`mask-before`)
    this._masterChart.xAxis[0].addPlotBand({
      id: `mask-before`,
      from: min,
      to: max,
      color: `rgba(0, 0, 0, 0.2)`
    })
    self._detailChart.showLoading(`<div class="spinner-border text-success" role="status">
                        <span class = "sr-only"> Loading...</span></div>`)
    this._axiosInstance.get(url).then(
      response => {
        const data = response.data.newChartData
        const sumToCheck = response.data.sumToCheck
        self._refLength = response.data.refLength
        console.log(self._detailChart.series)
        self._detailChart.series[0].setData(data.sequenced, false, false, false)
        self._detailChart.series[1].setData(data.unblocked, false, false, false)
        self._masterChart.xAxis[0].removePlotBand(`mask-before`)
        self._detailChart.hideLoading()
        // self._coverageChart.masterChart.hideLoading()
        self._detailChart.redraw()
        // self._coverageChart.masterChart.redraw()

        this._oldSumToCheck = sumToCheck
      }).catch(
      error => {
        console.error(error)
      })
    this._detailChart.xAxis[0].setExtremes(min, max)
  }
}
