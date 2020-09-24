class LiveMinKnowController {
  constructor (flowcellId) {
    this._flowcellId = flowcellId
    this._axiosInstance = axios.create({
      headers: { 'X-CSRFToken': getCookie(`csrftoken`) }
    })
    this._lastMinionRunStatsId = 0
    // Todo loop below to create charts
    this._liveYieldChart = makeLiveChart(
      `chart-yield`,
      `Yield Over Time`,
      `cumulative yield`
    )
    this._livePoreStatesChart = makeAreaPlot(
      `pore-states`,
      `Pore States`,
      `Pore States`
    )
    this._liveHistogramChart = makeLiveHistogram(
      `live-histogram`,
      `Histogram of Read Lengths (Events)`,
      `test`
    )
    this._liveOccupancyChart = makeLiveChart(
      `live-occupancy`,
      `% Occupancy Over Time`,
      `% Occupancy`
    )
    this._liveInStrandChart = makeLiveChart(
      `live-strand`,
      `Strands in pore counts`,
      `Number of pores in strand/single`
    )
    this._liveTemperatureChart = makeLiveChart(
      `live-temperature`,
      `Temperature Over Time`,
      `°Celcius`
    )

    this._liveVoltageChart = makeLiveChart(
      `live-voltage`,
      `Voltage Over Time`,
      `mV`
    )
    this._histogramHistory = null
    this._first = true
    if (getSelectedTab() === `live-event-data`) { this._makePageUnscrollable() }
    // this._interval = setInterval(this.fetchLiveEventsData, 60000, this._flowcellId, this)
    // $(window).on(`unload`, function () {
    //   console.log(`clearing base-called interval`)
    //   clearInterval(this._interval)
    // })
    $(`#loader-wrapper`).css(`width`, $(`#tab-live-event-data`).css(`width`))
  }

  /**
   * Prepare the histogram slider to change the displayed data by the histogram
   * @private
   */
  _prepareHistogramSlider () {
    const that = this
    $(`#histogram-date-picker`).on(`change`, (event) => {
      const index = parseInt($(`#${event.currentTarget.id}`).val())
      that._updateLiveHistogram(that._histogramHistory, index)
    })
  }

  /**
   * Make the page unscrollable whilst we are drawing the charts
   */
  _makePageUnscrollable () {
    $(`html`).addClass(`disable-scroll`)
  }

  /**
   * Reveal the live data tab by removing loder divs once charts are drawn
   */
  _revealPage () {
    $(`#tab-live-event-data`).addClass(`loaded`)
    $(`html`).removeClass(`disable-scroll`)
  }

  updateTab () {
    if (!$(`#tab-live-event-data`).hasClass(`loaded`)){
      this._makePageUnscrollable()
    }
    this.fetchLiveEventsData(this._flowcellId, this)
  }

  /**
   * Live histogram chart updated by us
   * @param histogramData {[{}]} list of series to display, starting with earliest sample time and values and
   * getting older, with each element
   * @param index {number} Chosen sample time to display
   *
   */
  _updateLiveHistogram (histogramData, index) {
    // find the bin index that has the n50 in it
    const n50 = histogramData[index].histogram_series.findIndex(e => e.n50 === true)
    $(`#histogram-date-picker`).attr({ max: `${histogramData.length - 1}`, value: `${index}` })
    $(`#displayed-date`).html(`${new Date(histogramData[index].sample_time).toGMTString()}`)
    this._liveHistogramChart.series[0].setData(histogramData[index].histogram_series)
    this._liveHistogramChart.xAxis[0].setCategories(histogramData[index].categories)
    this._liveHistogramChart.xAxis[0].removePlotBand(`plot-band-1`)
    this._liveHistogramChart.xAxis[0].addPlotBand({
      from: n50 - 0.5,
      to: n50 + 0.5,
      color: `#FCFFC5`,
      id: `plot-band-1`
    })
    this._liveHistogramChart.xAxis[0].removePlotBand(`plot-band-2`)
    this._liveHistogramChart.xAxis[0].addPlotBand({
      color: `black`,
      width: 2,
      dashStyle: `longdashdot`,
      value: n50,
      label: {
        text: `Estimated Read N50`,
        align: `left`,
        rotation: 0,
        x: +16 // Amount of pixels the label will be repositioned according to the alignment.
      },
      id: `plot-band-2`
    })
    // this._liveHistogramChart.xAxis[0].removePlotBand(`plot-band-3`)
    // this._liveHistogramChart.xAxis[0].addPlotBand({
    //   color: `black`,
    //   width: 2,
    //   dashStyle: `longdashdot`,
    //   value: (Math.floor(this.totalyield / this.readcount / this.datain2)),
    //   label: {
    //     text: `Estimated Read Average - ` + Math.round(this.totalyield / this.readcount / 1000 * 100) / 100 + ` K events`,
    //     align: `left`,
    //     rotation: 0,
    //     x: +10,
    //     y: +30 // Amount of pixels the label will be repositioned according to the alignment.
    //   },
    //   id: `plot-band-3`
    // })
    this._liveHistogramChart.reflow()
  }

  /**
   * Update given chart on the live tab with the provided data. Checks if the data is different, and replaces
   * already extant data.
   * @param data {[{}]} a list of series objects to be inserted
   * @param chart {} Highcharts Api object
   * @param runStarts [[]] Array of Arrays like [[run_id, run_start_time]], where run_start_time is in microseconds
   */
  updateLiveTabChart (data, chart, runStarts) {
    let redraw = false
    data.forEach(chartReadySeries => {
      console.log(chartReadySeries)
      // check if this series has already been added
      if (chart.series.findIndex(e => e.name === chartReadySeries.name) + 1) {
        const seriesToUpdate = chart.series.find(e => e.name === chart.name)
        // check if the new series has a larger length
        if (chartReadySeries.length > seriesToUpdate.length) {
          const differenceInLength = chartReadySeries.length - seriesToUpdate.length
          if (differenceInLength < 150) {
            // if we are only adding a few new points, let's add them into existing series
            chartReadySeries.data.forEach(point => {
              seriesToUpdate.addPoint(point, false)
            })
            redraw = true
          } else {
            seriesToUpdate.setData(chartReadySeries.data, false)
            redraw = true
          }
        }
      } else {
        console.log(`addIngSeties`)
        chart.addSeries(chartReadySeries, false)
        // runStarts.forEach(([runName, runStart]) => {
        //   chart.xAxis[0].addPlotLine({
        //     label: {
        //       text: `Run: ${runName}`,
        //       x: 10,
        //       verticalAlign: `top`,
        //       style: {
        //         color: `blue`,
        //         fontWeight: `bold`,
        //         fontSize: `.5rem`
        //       }
        //     },
        //     color: `red`,
        //     value: runStart,
        //     dashStyle: `longdashdot`,
        //     width: 2
        //   })
        // })
        redraw = true
      }
    })
    if (redraw) { console.log(chart); chart.redraw() }
  }

  fetchLiveEventsData (flowcellId, that) {
    if (getSelectedTab() !== `live-event-data`) {
      return
    }
    that._axiosInstance.get(`/api/v1/reads/flowcells/${flowcellId}/runstats/${that._lastMinionRunStatsId}`).then(
      response => {
        if (response.status === 206) { return }
        const index = that._first ? response.data.histogram_history.length - 1 : parseInt($(`#histogram-date-picker`).val())
        const runInfo = response.data.run_info
        that._first = false
        that._histogramHistory = response.data.histogram_history
        that.updateLiveTabChart(response.data.yield_history, that._liveYieldChart, runInfo)
        console.log(response.data.pore_history)
        that.updateLiveTabChart(response.data.pore_history, that._livePoreStatesChart, runInfo)
        that._prepareHistogramSlider()
        that._updateLiveHistogram(response.data.histogram_history, index)
        that.updateLiveTabChart(response.data.occupancy_history, that._liveOccupancyChart, runInfo)
        that.updateLiveTabChart(response.data.in_strand_history, that._liveInStrandChart, runInfo)
        that.updateLiveTabChart(response.data.temperature_history, that._liveTemperatureChart, runInfo)
        that.updateLiveTabChart(response.data.voltage_history, that._liveVoltageChart, runInfo)
        that._revealPage()
        that._lastMinionRunStatsId = response.data.last_minion_run_stats_id
      }
    )
  }
}
