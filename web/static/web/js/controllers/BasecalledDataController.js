/**
 * @author Adonis
 */
class BasecalledDataController {
    /**
     *
     * @param flowcellId {string} The PK of the flowcell record in the database.
     */
    constructor(flowcellId) {
        // TODO BUG could clash if multiple browser tabs open
        setSelectedBarcode('All reads', "BasecalledData");
        this._currentBarcode = getSelectedBarcode("BasecalledData");
        this._flowcellId = flowcellId;
        this._axiosInstance = axios.create({
            headers: {"X-CSRFToken": getCookie('csrftoken')}
        });
        this._barcodesHtmlListElement = $("#nav-tabs-barcodes");
        this.barcodesList = [];
        this.barcodesElementsList = [];
        this._basecalledSummaryTable = null;
        this._updateBarcodeNavTab(flowcellId);
        this._fetchSummaryDataHtmlTable(flowcellId);
        this._createHistogramCharts(flowcellId, this._currentBarcode);
        this._createBaseCalledReadCharts(flowcellId, this._currentBarcode);
        this._createFlowcellHeatMaps(flowcellId);
        this._interval = setInterval(this._updateTab, 60000, flowcellId, this._currentBarcode, this, false);
        $(window).on("unload", function () {
            console.log("lcearing interval");
            clearInterval(this._interval);
        });
    }

    /**
     * Change class on the loading div to reveal results.
     * @private
     */
    _revealResults(){
        $("#tab-basecalled-data").addClass("loaded");
    }

    /**
     * Create and add any new barcode tabs to the barcode nav tab.
     * @param flowcellId {string} The primary key of the flowcell record in the database.
     * @private
     */
    _updateBarcodeNavTab(flowcellId) {
        let firstIteration, that, index, barcodes,
            barcode_objs, active, newBarcodeElementsArray = [];
        that = this;
        firstIteration = this.barcodesList.length === 0;
        this._axiosInstance.get("/api/v1/flowcells/barcodes", {params: {flowcellId}}).then((response) => {
            if (response.status !== 200) console.error(`Error ${response.status}`);
            barcode_objs = response.data;
            // Remove no barcode
            barcode_objs = barcode_objs.filter(el => {
                return el.name !== "No barcode";
            })
            // Just string names as an array
            barcodes = barcode_objs.map(barc => barc.name)
            // If we already have HTML elements
            if (!firstIteration) {
                that.barcodesList.forEach(presBarcode => {
                    // if the barcode in controller in the list of barcodes fetched from the server
                    if (barcodes.includes(presBarcode)) {
                        // get the index,
                        index = barcode_objs.findIndex(b => b.name === presBarcode);
                        // remove it from the freshly fetched barcodes
                        barcode_objs.splice(index, 1);
                    }
                });
            }
            barcode_objs.forEach(barcode => {
                that.barcodesList.push(barcode.name)
                active = barcode.name === that._currentBarcode ? "active" : "";
                // that.barcodesElementsList.push(`<li class="barcode-tab nav-item">
                //                             <a class="nav-link ${active}" data-toggle="pill" id="${barcode.name.replace(" ", "_")}" onclick="basecalledDataController._changeBarcode(event)">${barcode.name}</a>
                //                         </li>`)
                newBarcodeElementsArray.push(`<li class="barcode-tab nav-item">
                                            <a class="nav-link ${active}" data-toggle="pill" id="${barcode.name.replace(" ", "_")}" onclick="basecalledDataController._changeBarcode(event)">${barcode.name}</a>
                                        </li>`)
            })
            // this.barcodesElementsList = this.barcodesElementsList.sort()
            if (barcode_objs.length > 0) {
                this._barcodesHtmlListElement.append(newBarcodeElementsArray)
            }
        }).catch(error => {
            console.error(error)
        });
    };

    /**
     * Change the active barcode. Add active class to new barcode, remove from old, call function to draw new charts.
     * @param event - Event object from event listener
     * @private
     */
    _changeBarcode(event) {
        console.log("change BArcode")
        let selectedBarcode = event.target.innerText;
        // Remove active from the selected barcode class
        $(`#${this._currentBarcode.replace(" ", "_")}`).removeClass("active")
        setSelectedBarcode(selectedBarcode, "BasecalledData")
        this._currentBarcode = selectedBarcode;
        event.target.classList.add("active")

        // call the function to draw all the charts here.
        this._updateTab(this._flowcellId, selectedBarcode, this, true)
    }

    /**
     * Request the HTML summary tables to be appended to the top of the Basecalled data tabs.
     * @param flowcellId {number} Primary key of the flowcell record in the mysql database.
     * @private
     */
    _fetchSummaryDataHtmlTable(flowcellId) {
        let that = this;
        this._basecalledSummaryTable = $('#basecalled-summary');

        this._axiosInstance.get(`/flowcells/${flowcellId}/basecalled-summary-html`).then((response) => {
            if (response.status !== 200) {
                console.error(`Error, incorrect status, expected 200, got ${response.status}`)
            }
            that._basecalledSummaryTable.html(response.data)
            setTimeout(this._revealResults, 100);

        }).catch(error => {
            console.error(error)
        })
    }

    /**
     * Check if a given chart is already initialised, and if not, initialise it
     * @param chart {obj} Chart object
     * @param divId {string} Id of the div, to initialise chart in.
     * @param chartTitle {string} Title of the chart.
     * @param yAxisTitle {string} Title of the y Axis.
     * @param collect {bool} Whether chart is proportional or not.
     * @param categories {array} List of categories for histogram.
     * @param heatmap {boolean} If chart to be checked is of type heatmap.
     * @private
     */
    _checkChartExists(chart, divId, chartTitle, yAxisTitle, collect, categories, heatmap = false) {
        // highcharts options
        let options = {};
        if (categories.length > 0) {
            options = {
                plotOptions: {
                    column: {
                        stacking: 'normal'
                    }
                },
                chart: {
                    type: 'column'
                },
                xAxis: {
                    type: 'category',
                    categories: categories
                }
            };
            if (collect) {
                options["yAxis"] = {"max": 100}
            }
        }

        // check if the chart is referenced by the class already. Should be if chart is initialised.
        if (!this.hasOwnProperty(chart)) {
            if (!heatmap) {
                chart = makeSplineChart(divId, chartTitle, yAxisTitle)
                chart.update(options)
            } else if (heatmap) {
                chart = makeHeatmapChart(divId, chartTitle)
            }
        }
        return chart
    }

    /**
     * Create highChart elements for the four histogram charts on the base-called data tab.
     * @param flowcellId {string} The primary key of the flowcell record in the database.
     * @param barcodeName {string} The name of the selected barcode.
     * @private
     */
    _createHistogramCharts(flowcellId, barcodeName) {
        /*
         * Request histogram data - Request data for the two histograms on the basecalled data tab
         */
        let initialise = true;
        this._axiosInstance.get(`/api/v1/flowcells/${flowcellId}/histogram-summary/`, {
            params: {
                barcodeName,
                initialise
            }
        }).then(response => {
            let categories = response.data;
            let charts = [{
                "chart": "_chartHistogramBasesSequencedByReadLength",
                "divId": "histogram-bases-sequenced-by-read-length",
                "chartTitle": "HISTOGRAM OF BASES SEQUENCED BY READ LENGTH",
                "yAxisTitle": "NUMBER OF BASES",
                "collect": false
            }, {
                "chart": "_chartHistogramReadLength",
                "divId": "histogram-read-lengths",
                "chartTitle": "HISTOGRAM OF READ LENGTHS",
                "yAxisTitle": "NUMBER OF READS",
                "collect": false
            }, {
                "chart": "_collectchartHistogramBasesSequencedByReadLength",
                "divId": "collect-histogram-bases-sequenced-by-read-length",
                "chartTitle": "COLLECT HISTOGRAM OF BASES SEQUENCED BY READ LENGTH",
                "yAxisTitle": "PROPORTION OF BASES LONGER THAN X",
                "collect": true
            }, {
                "chart": "_collectchartHistogramReadLength",
                "divId": "collect-histogram-read-lengths",
                "chartTitle": "COLLECT HISTOGRAM OF READ LENGTHS",
                "yAxisTitle": "PROPORTION OF READS LONGER THAN X",
                "collect": true
            }]
            if (response.status !== 200) {
                console.error(`Error, incorrect status, expected 200, got ${response.status}`)
            }

            charts.forEach(chart => {
                // set the this chart pointer value of the class to the highcharts object returned
                this[chart.chart] = this._checkChartExists(...Object.values(chart), categories)


            });
            this._updateHistogramChartsData(flowcellId, barcodeName)
        }).catch(error => {
            console.error(error)
        })
    }

    /**
     * Update highChart data for the four histogram charts on the base-called data tab.
     * @param flowcellId {string} The primary key of the flowcell record in the database.
     * @param barcodeName {string} The name of the selected barcode.
     * @param changeBarcode {boolean} If update is triggered by barcode change.
     * @private
     */
    _updateHistogramChartsData(flowcellId, barcodeName, changeBarcode) {
        let charts = [["read_length", this._chartHistogramBasesSequencedByReadLength],
            ["read_count", this._chartHistogramReadLength],
            ["collect_read_length", this._collectchartHistogramBasesSequencedByReadLength],
            ["collect_read_count", this._collectchartHistogramReadLength]]
        let preExistingSeriesNames;
        let chartData;

        this._axiosInstance.get(`/api/v1/flowcells/${flowcellId}/histogram-summary/`, {params: {barcodeName}}).then(response => {
            let chartPoints = response.data.data;
            if (response.status !== 200) {
                console.error(`Error, incorrect status, expected 200, got ${response.status}`)
            }
            charts.forEach(chart => {
                let newDataKey = chart[0]
                let chartReference = chart[1]
                let difference; // Difference between the names of the newly fetched data and the data in the charts
                let union; // new data for a series we have already drawn
                if (changeBarcode) {
                    clearChartData(chartReference)
                }
                preExistingSeriesNames = chartReference.series.map(e => e.name)

                chartData = chartPoints[newDataKey];
                difference = chartData.filter(x => !preExistingSeriesNames.includes(x.name));
                // add any new series
                difference.forEach(newSeries => {
                    chartReference.addSeries(newSeries)
                })

                union = chartData.filter(x => preExistingSeriesNames.includes(x.name));
                union.forEach(updateToSeries => {
                    let seriesToUpdate = chartReference.series.filter(e => {
                        return e.name === updateToSeries.name
                    })
                    if (checkHighChartsDataIsNew(updateToSeries.data, seriesToUpdate[0].options.data)) {
                        console.log("data is identical skipping redraw.")
                    } else {
                        seriesToUpdate[0].setData(updateToSeries.data)
                    }
                })

            })
        }).catch(error => {
            console.error(error)
        })
    }

    /**
     * Create highchart objects for the six line charts on the base-called data tab.
     * @param flowcellId {string} The primary key of the flowcell record in the database.
     * @param barcodeName {string} The name of the selected barcode.
     * @private
     */
    _createBaseCalledReadCharts(flowcellId, barcodeName) {


        let charts = [{
            "chart": "_averageQualityOverTime",
            "divId": "average-quality-over-time",
            "chartTitle": "AVERAGE QUALITY OVER TIME",
            "yAxisTitle": "AVERAGE READ QUALITY SCORE",
        }, {
            "chart": "_averageReadLengthsOverTime",
            "divId": "average-read-lengths-over-time",
            "chartTitle": "AVERAGE READ LENGTH OVER TIME",
            "yAxisTitle": "AVERAGE READ LENGTH",
        }, {
            "chart": "_cumulativeYieldOverTime",
            "divId": "cumulative-yield-over-time",
            "chartTitle": "CUMULATIVE BASES",
            "yAxisTitle": "CUMULATIVE BASES",
        }, {
            "chart": "_cumulativeNumberReadsOverTime",
            "divId": "cumulative-number-reads-over-time",
            "chartTitle": "CUMULATIVE READS",
            "yAxisTitle": "CUMULATIVE READS",
        }, {
            "chart": "_maxReadLengthsOverTime",
            "divId": "max-read-lengths-over-time",
            "chartTitle": "MAX READ LENGTH OVER TIME",
            "yAxisTitle": "MAX READ LENGTH",
        }, {
            "chart": "_sequencingRate",
            "divId": "sequencing-rate",
            "chartTitle": "SEQUENCING RATE",
            "yAxisTitle": "BASES/SECOND",
        }]

        charts.forEach(chart => {
            // set the this chart pointer value of the class to the highcharts object returned
            this[chart.chart] = this._checkChartExists(...Object.values(chart), false, [])
        })
        this._updateBaseCalledReadCharts(flowcellId, barcodeName)
    }

    //TODO add in the run start lines.
    /**
     * update highchart data for the six line charts on the base-called data tab.
     * @param flowcellId {string} The primary key of the flowcell record in the database.
     * @param barcodeName {string} The name of the selected barcode.
     * @param changeBarcode {boolean} If update is triggered by barcode change.
     * @private
     */
    _updateBaseCalledReadCharts(flowcellId, barcodeName, changeBarcode) {
        let charts = [this._averageQualityOverTime,
            this._averageReadLengthsOverTime,
            this._cumulativeYieldOverTime,
            this._cumulativeNumberReadsOverTime,
            this._maxReadLengthsOverTime,
            this._sequencingRate]
        console.log("Updating line charts");
        this._axiosInstance.get(`/api/v1/flowcells/${flowcellId}/statistics/`, {params: {barcodeName}}).then(response => {
            let key, value;
            //TODO add in the run start lines.
            let runData = response.data["runs"];
            let chartData = response.data.data;
            let newChartDataArray;
            let newChartData;
            let chartCleared = false; // clear the chart once for old data on first new series addition
            let identical = false; // If we have no new data
            let chartHasData = false; // If chart has no data
            let seriesToUpdate = {}; // The series in the chart we are updating
            console.time("statistics charts")
            if (response.status !== 200) {
                console.error(`Error, incorrect status, expected 200, got ${response.status}`)
            }
            // for a series worth of data
            for ([key, value] of Object.entries(chartData)) {
                // map together our super weird data structure, cry
                newChartDataArray = [value.map(x => [x[0], x[1]]),
                    value.map(x => [x[0], x[2]]),
                    value.map(x => [x[0], x[3]]),
                    value.map(x => [x[0], x[4]]),
                    value.map(x => [x[0], x[5]]),
                    value.map(x => [x[0], x[6]])]

                charts.forEach((chart, index) => {
                    // TODO this is where we would switch if we wanted multiple barcodes on a run
                    if (changeBarcode && !chartCleared) {
                        clearChartData(chart)
                    }
                    identical = false;
                    chartHasData = chart.series.length > 0
                    newChartData = newChartDataArray[index]
                    // See if this series is already in the chart
                    seriesToUpdate = chart.series.filter(x => {
                        return x.name === key
                    });
                    // only check if there is already data in the same chart series
                    if (chartHasData && seriesToUpdate.length > 0) {
                        identical = checkHighChartsDataIsNew(newChartData, seriesToUpdate[0].options.data)
                    }
                    if (identical) {
                        console.log("Data is the same. Skipping update.")
                        return
                    } else {
                        // only update if there is already data in the same chart series
                        if (chartHasData && seriesToUpdate.length > 0) {
                            seriesToUpdate[0].setData({
                                name: key,
                                data: newChartData
                            })
                        } else {
                            chart.addSeries({
                                name: key,
                                data: newChartData
                            })
                        }
                    }
                })
                chartCleared = true;
            }
            console.timeEnd("statistics charts")
        }).catch(error => {
            console.error(error)
        });
    }

    /**
     * Create the flowcell heatmaps that are at the bottom of the basecalled data page.
     * @param flowcellId {string} The primary key of the
     * @private
     */
    _createFlowcellHeatMaps(flowcellId) {
        let charts = [{
            "chart": "_chartReadsPerPore",
            "divId": "reads-per-pore",
            "chartTitle": "READS PER CHANNEL",
        }, {
            "chart": "_chartBasesPerPore",
            "divId": "bases-per-pore",
            "chartTitle": "BASES (KB) PER CHANNEL",
        }]

        charts.forEach(chart => {
            // set the this chart scope value of the class to the highcharts object returned
            this[chart.chart] = this._checkChartExists(...Object.values(chart), "", false, [], true)
        })

        this._updateFlowcellHeatMaps(flowcellId)
    }

    /**
     * Update the flowcell heatmaps on the basecalled data tab.
     * @param flowcellId {string} The primary key of the flowcell record in the mysql database.
     * @private
     */
    _updateFlowcellHeatMaps(flowcellId) {
        let charts = [["readCount", this._chartReadsPerPore], ["poreYield", this._chartBasesPerPore]];
        let chartReadyData = {"readCount": [], "poreYield": []}
        let chartData;
        let chartHasData;
        let size; // size of the flowcell. 512 for minIon etc.
        this._axiosInstance(`/api/v1/flowcells/${flowcellId}/channel-summary`).then(response => {
            console.log(response.data)
            if (response.status !== 200) {
                console.error(`Error, incorrect status, expected 200, got ${response.status}`)
            }
            chartData = response.data;
            size = response.data.length
            // again manipulate our data into a different format, cry
            chartData.forEach(channelDataPoint => {
                chartReadyData["readCount"].push([channelDataPoint[0], channelDataPoint[1], channelDataPoint[2]])
                chartReadyData["poreYield"].push([channelDataPoint[0], channelDataPoint[1], channelDataPoint[3]])
            })
            // add data to the charts!
            charts.forEach(([key, chart]) => {
                chartHasData = chart.series.length > 0;
                if (chartHasData) {
                    if (checkHighChartsDataIsNew(chartReadyData[key], chart.series[0].options.data)) {
                        console.log("Data is identical, skipping redraw.");
                    } else {
                        chart.series[0].setData(chartReadyData[key])
                    }
                } else {
                    chart.addSeries({
                        name: key,
                        data: chartReadyData[key],
                        borderWidth: size === 512 ? 1 : 0.5,
                    });
                }
            })
        }).catch(error => {
            console.error(error)
        });
    }

    /**
     * Update the data, called in interval in initialisation.
     * @param flowcellId {string} Primary key of the flowcell record in the mysql database.
     * @param barcodeName {string} Chose barcode of the tab.
     * @param that {BasecalledDataController} The scope of the class
     * @param changeBarcode {boolean} If update is triggered by barcode change.
     * @private
     */
    _updateTab(flowcellId, barcodeName, that, changeBarcode) {
        console.log("Updating!");
        console.log(barcodeName)
        let barcode = getSelectedBarcode("BasecalledData")
        that._fetchSummaryDataHtmlTable(flowcellId)
        that._updateHistogramChartsData(flowcellId, barcode, changeBarcode);
        that._updateBaseCalledReadCharts(flowcellId, barcode, changeBarcode);
        that._updateBarcodeNavTab(flowcellId, changeBarcode);
    }
}