// function updateLiveCumuYield ( liveTemperature, liveVoltage, livedata) {
//   // if (liveCumuYield.series.length < 1) {
//   //   liveCumuYield.addSeries({
//   //     data: livedata.yield_history
//   //   })
//   //   liveCumuYield.series[0].update({ name: `Events` }, false)
//   // } else {
//   //   liveCumuYield.series[0].setData(livedata.yield_history)
//   //   liveCumuYield.series[0].update({ name: `Events` }, false)
//   // }
//   //
//   // liveCumuYield.redraw()
//   // liveCumuYield.reflow()
//
//
//   // liveOccupancy.redraw()
//
//
//
//
//
//
//   liveVoltage.series[0].update({ name: `Voltage` }, false)
//
//   liveVoltage.redraw()
//   liveVoltage.reflow()
// }
//
// /**
//  * Make the page unscrollable whilst we are drawing the charts
//  */
// function makePageUnscrollable () {
//   $(`html`).addClass(`disable-scroll`)
// }
// /**
//  * Reveal the live data tab by removing loder divs once charts are drawn
//  */
// function revealPage () {
//   $(`#tab-live-event-data`).addClass(`loaded`)
//   $(`html`).removeClass(`disable-scroll`)
// }
//
// function requestLiveRunStats (id) {
//   console.log(`This is lastRead from flowcellcontroller: ` + flowcellController.lastRead)
//
//   // var url_livestats = '/api/v1/flowcells/' + id + '/runstats/' + flowcell_controller.lastRead;
//   var url_livestats = `/api/v1/flowcells/` + id + `/runstats/0`
//
//   $.get(url_livestats, function (result) {
//     // var liveHistogram = makeLiveHistogram(
//     //   `live-histogram`,
//     //   `Histogram of Read Lengths (Events)`,
//     //   `test`
//     // )
//
//     // var liveCumuYield = makeLiveChart(
//     //   `chart-yield`,
//     //   `Yield Over Time`,
//     //   `cumulative yield`
//     // )
//
//
//
//
//
//     var data = result.data
//     //
//     var live_data = new Array()
//     live_data.colours_string = result.minKNOW_colours_string
//     live_data.voltage = new Array()
//     live_data.asictemp = new Array()
//     live_data.heatsinktemp = new Array()
//     live_data.strand = new Array()
//     live_data.adapter = new Array()
//     live_data.good_single = new Array()
//     live_data.pore = new Array()
//     live_data.currpercentage = null
//     live_data.currstrand = null
//     live_data.percentage = new Array()
//     live_data.yield_history = new Array()
//     live_data.meanratio_history = new Array()
//     live_data.instrand_history = new Array()
//     live_data.openpore_history = new Array()
//     live_data.pore_history = new Array()
//
//     // var myStringArray = [`above`, `adapter`, `below`, `good_single`, `strand`, `inrange`, `multiple`, `pending_mux_change`, `saturated`, `unavailable`, `unblocking`, `unclassified`, `unknown`, `zero`, `pore`, `no_pore`]
//     //
//     // var arrayLength = myStringArray.length
//     //
//     // for (var i = 0; i < arrayLength; i++) {
//     //   live_data.pore_history[myStringArray[i]] = new Array()
//     // }
//
//     if (data.length > 0) {
//       data.sort(function (a, b) {
//         var a_date = new Date(a.sample_time)
//         var b_date = new Date(b.sample_time)
//
//         return a_date - b_date
//       })
//
//       for (var i = 0; i < data.length; i++) {
//         timestamp = new Date(data[i].sample_time).getTime()
//         live_data.live_read_count = data[i].minKNOW_read_count
//         live_data.voltage.push([timestamp, data[i].voltage_value])
//         live_data.asictemp.push([timestamp, data[i].asic_temp])
//         live_data.heatsinktemp.push([timestamp, data[i].heat_sink_temp])
//         live_data.strand.push([timestamp, data[i].strand])
//         live_data.adapter.push([timestamp, data[i].adapter])
//         live_data.pore.push([timestamp, data[i].pore])
//         live_data.good_single.push([timestamp, data[i].good_single])
//         live_data.currpercentage = data[i].occupancy
//         live_data.currstrand = data[i].strand
//         live_data.percentage.push([timestamp, data[i].occupancy])
//         // live_data.yield_history.push([timestamp, data[i].event_yield])
//         live_data.meanratio_history.push([timestamp, data[i].mean_ratio])
//         live_data.instrand_history.push([timestamp, data[i].in_strand])
//         live_data.openpore_history.push([timestamp, parseInt(data[i].open_pore)])
//       }
//
//       // calculatereadtoeventscaling();
//
//       // updateLiveHistogram(liveHistogram, data)
//       // updateLiveCumuYield(  liveTemperature, liveVoltage, live_data)
//       revealPage()
//     }
//   })
// };
