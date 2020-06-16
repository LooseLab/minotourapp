/**
 * TODO this would be a good candidate to list the run the data was coming from in localStorage
 */
class ArticController {
    constructor(flowcellId) {
        this._articCoverageScatterMaster = "hello";
        this._articCoverageScatterDetail = "world";
        this._flowcellId = flowcellId;
        this._first = true;
        this._showWidth = 250000;
        this._barcodeChosen = null;
        this._logged = $("#log-coverage")[0].checked ? 1 : 0;
        this._axiosInstance = axios.create({
            headers: {"X-CSRFToken": getCookie('csrftoken')}
        });
        this._perBarcodeAverageReadLengthChart = null;
        this._perBarcodeCoverageChart = null;
        this._coverageSummaryTable = $("#artic-coverage-summary-table");
        this._articProportionPieChart = null;
        this._drawArticCharts(this);
        this._interval = setInterval(this._drawArticCharts, 45000, this);
        $(window).on("unload", function () {
            console.log("clearing Artic interval");
            clearInterval(this._interval);
        });

    }

    /**
     *  Update the series of an existing chart. Checks to see if data is identical
     * @param {obj} chart - Chart object that is being updated
     * @param {Array.Array} newSeries - The new data to add into the chart
     * @param {number} index - The index of the series that we are updating. Use 0 .
     */
    _updateExistingChart(chart, newSeries, index) {
        // Set the data on an existing series on a chart, either master or detail
        // if the data is identical
        if (checkHighChartsDataIsNew(newSeries, chart.series[index].options.data)) {
            console.log("Data is identical , skip redraw");
        } else {
            chart.showLoading('Fetching data from the server');
            chart.series[index].setData(newSeries);
            chart.hideLoading();
        }
    }

    /**
     * Render the smaller coverage master chart that shows the region selected, and allows users to select a region.
     * @param data {object} Object conatinf response data from the server.
     * @param flowcellId {number} The Primary key of the flowcell record in the database.
     */
    _createCoverageMaster(data, flowcellId) {
        const self = this;
        let that = this;
        this._axiosInstance.get("/api/v1/artic/fetch-amplicon-bands").then(response => {
            let bands = response.data.amplicon_band_coords
            let colourLookup = response.data.colours
            let highChartsBands = []
            let detailsBands = []
            bands.forEach((band, index) => {
                let masterBand = {
                    from: band[0],
                    to: band[1],
                    color: colourLookup[band[2]],
                }
                highChartsBands.push(masterBand)
                detailsBands.push({
                    ...masterBand, ...{
                        label: {
                            text: `# ${index + 1}`,
                            style: {
                                color: "grey"
                            },
                            rotation: 270
                        }
                    }
                })
            })

            this._articCoverageScatterMaster = Highcharts.chart("coverageArticMaster", {
                chart: {
                    zoomType: "x",
                    height: "150px",
                    marginLeft: 50,
                    marginRight: 20,
                    events: {
                        selection: that._afterSelection.bind(this)
                    }

                },
                boost: {
                    useGPUTranslations: true,
                    usePreAllocated: true
                },
                xAxis: {
                    min: 0,
                    max: data.coverage.xmax,
                    gridLineWidth: 1,
                    zoomEnabled: false,
                    maxZoom: 10000,
                    title: {
                        text: "Base Number"
                    },
                    plotBands: [...highChartsBands, {
                        id: 'mask-before',
                        from: data.coverage.data[0][0],
                        to: data.coverage.data[data.coverage.data.length - 1][0],
                        color: 'rgba(0, 0, 0, 0.2)'
                    }]
                },
                yAxis: {
                    min: 0,
                    max: data.coverage.ymax,
                    title: {
                        text: "Coverage"
                    }
                },
                legend: {
                    enabled: false
                },
                title: {text: null},
                plotOptions: {
                    series: {
                        fillColor: {
                            linearGradient: [0, 0, 0, 70],
                            stops: [
                                [0, Highcharts.getOptions().colors[0]],
                                [1, 'rgba(255,255,255,0)']
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
                    type: "area",
                    color: "blue",
                    data: data.coverage.data,
                    marker: {
                        fillcolor: "red",
                        radius: 1
                    },
                    tooltip: {
                        followPointer: false,
                        pointFormat: "[{point.x:.1f}, {point.y:.1f}]",
                        valueDecimals: 8
                    }
                }]
                // Call back function, called after the master chart is drawn, draws the detail charts
            }, () => {
                console.log(detailsBands)
                that._articCoverageScatterDetail = that._createDetailChart("coverageArticDetail", data.coverage.ymax, "Coverage", "line"
                    , false, false, detailsBands);
            });
        }).catch(error => {
            console.error(error)
        })
    }

    /**
     * Render the coverage detail chart, with empty series
     * @param {string} targetDivId Div ID chart is drawn to
     * @param {number} ymax The maximum value of the Y axis
     * @param {string} title The chart title
     * @param {string} chartType Chart title
     * @param {boolean} stepped Stepped line
     * @param {boolean} legend Display chart legend
     * @param {Array} bands Array of plot bands for amplicons
     * @returns {obj} chartyChart Highcharts chart.
     */
    _createDetailChart(targetDivId, ymax, title, chartType = "scatter", stepped = false, legend = false, bands = []) {
        // The data we are going to use to display, starts a display with an empty chart
        let detailData = [];
        // The variable assigned to the chart we will return
        let chartyChart;

        let that = this;

        chartyChart = Highcharts.chart(targetDivId, {

            chart: {
                zoomType: "x",
                height: "300px",
                reflow: true,
                events: {
                    selection: that._afterSelection.bind(this)
                }
            },

            // boost: {
            //     useGPUTranslations: true,
            //     usePreAllocated: true
            // },

            xAxis: {
                type: chartType,
                gridLineWidth: 1,
                title: {
                    text: "Base number"
                },
                plotBands: bands
            },

            yAxis: {
                // Renders faster when we don"t have to compute min and max
                // type: 'logarithmic',
                softMin: 0,
                max: ymax,
                title: {
                    text: null
                },
                startOnTick: false,
                endOnTick: true
            },

            title: {
                text: title
            },

            legend: {
                enabled: legend
            },
            tooltip: {
                formatter: function () {
                    return `The value for base <b>${this.x}</b> is <b>${this.y}</b>`;
                }
            },
            exporting: {"enabled": true},

            series: [{
                type: chartType,
                color: "blue",
                step: stepped,
                data: detailData,
                marker: {
                    fillcolor: "red",
                    radius: 1
                },
                tooltip: {
                    followPointer: false,
                    pointFormat: "[{point.x:.1f}, {point.y:.1f}]",
                    valueDecimals: 8
                }
            }]

        });
        return chartyChart;
    }

    /**
     * The event that is called after an area selection is made on the master or detail chart. Sets the plot band, calls to server for new data.
     * @param event {object} The event object generated when the event is triggered
     */
    _afterSelection(event) {
        // Event function called after the zoom box is chosen
        let label;
        let that = this;
        let flowcellId = this._flowcellId;
        let charts = [this._articCoverageScatterDetail];
        // The minimum and maximum value on the X axis
        let min = Math.trunc(event.xAxis[0].min);
        let max = Math.trunc(event.xAxis[0].max);
        let logCoverage;
        // Prevent the default behaviour
        event.preventDefault();
        //
        // Change the grey plot band on the master chart to reflect the selected region
        that._articCoverageScatterMaster.xAxis[0].removePlotBand('mask-before');
        // Add the new one in the correct location
        that._articCoverageScatterMaster.xAxis[0].addPlotBand({
            id: 'mask-before',
            from: min,
            to: max,
            color: 'rgba(0, 0, 0, 0.2)'
        });

        // Show the loading screen on the charts
        logCoverage = this._logged ? 1 : 0;
        charts.forEach(function (chart) {
            chart.showLoading('Fetching data from the server');
        });
        this._axiosInstance.get('/api/v1/artic/visualisation/detail', {
            params: {min, max, barcodeChosen: this._barcodeChosen, flowcellId, logCoverage}
        }).then(response => {
            console.log("response in selecetion")
            // if successful request, set the detail charts to show the newly retrieve data
            if (response.status === 200) {
                that._updateAxesAndData(that._articCoverageScatterDetail, min, max, response.data.coverage, true);
            } else if (response.status === 404) {
                console.error("Artic detail chart data not found.")
            } else {
                throw "Request error";
            }
        }).catch(error => {
            console.error(error)
        });
    }

    /**
     * Update the axes and series data of a specified chart.
     * @param chart {object} The highcharts APi object for the detail and the master chart.
     * @param min {number} The x axis minimum value.
     * @param max {number} The x axis maximum value.
     * @param data {Array[]} Array of Arrays, specifying cartesian coordinates for the plot.
     * @param ymax {number} The maximum value of the y axis
     * @private
     */
    _updateAxesAndData(chart, min, max, data, ymax) {
        let yAxisType;
        // Y max let's us know we need to change the yAxis
        // Set the xAxis extremes to the newly selected values
        chart.xAxis[0].setExtremes(min, max);
        // if we need to update the yaxis as well
        if (ymax) {
            // set the yAxisextremes, if it's logged set to 0.1
            if (this._logged) {
                yAxisType = "logarithmic"
                chart.yAxis[0].setExtremes(0.1, data.ymax)
                chart.update({
                    yAxis: {type: yAxisType}
                })
            } else {
                yAxisType = "linear"
                chart.update({
                    yAxis: {type: yAxisType}
                })
                chart.yAxis[0].setExtremes(0, data.ymax)
            }
            console.log("updating log type");
            console.log(yAxisType)
            console.log(chart)

        }
        this._data = data
        // set the data for the chart to display
        chart.series[0].setData(data.data);
        // update the axes
        chart.xAxis[0].update();

        // hide the loading screen
        chart.hideLoading();
    }

    /**
     * Draw the select box on the artic page for the coverage tracker charts.
     * @param flowcellId {number} Primary key of the flowcell that barcodes are bring added to the selection box for.
     * @private
     */
    _createOrUpdateSelectionBox(flowcellId) {
        let firstIteration, optionChildren, barcodes, alreadyPresentBarcode, index, that, $current, barcodeList;
        that = this;
        barcodeList = [];
        // for each select box on the page
        $current = $("#select-artic-barcode");
        optionChildren = $current.children('option');
        firstIteration = optionChildren.length === 0 ? true : false;
        that._axiosInstance.get("/api/v1/artic/visualisation/barcodes", {params: {flowcellId}}).then((response) => {
            if (response.status !== 200) console.error("Error");
            barcodes = response.data;
            // else we need to remove the barcodes that we have already added then add them.
            if (!firstIteration) {
                alreadyPresentBarcode = [...$("#select-artic-barcode")[0].children];
                alreadyPresentBarcode.forEach(presBarcode => {
                    // if the chromsome in already present chromsome is in the list of barcodes fetched from the server
                    if (barcodes.includes(presBarcode.textContent)) {
                        // get the index,
                        index = barcodes.indexOf(presBarcode.textContent);
                        // remove it from the freshly fetched barcodes
                        barcodes.splice(index, 1);
                    }
                });
            } else {
                barcodes.unshift("Choose Barcode");
            }
            barcodes.forEach(barcode => {
                if (barcode === "Choose Barcode") {
                    barcodeList.push(`<option class="sel__box__options" id="place-holder" value="${barcode}">${barcode}</option>`)
                } else {
                    barcodeList.push(`<option class="sel__box__options" value="${barcode}">${barcode}</option>`)
                }
            })
            if (barcodeList.length) {
                $current.html(barcodeList.join(""))
            }
            $current.change(function () {
                let txt;
                txt = $(this).val();
                console.log(txt)
                that._barcodeChosen = txt;
                // Remove the placeholder "Choose Barcode" option
                $("#place-holder").remove()
                // enable the log coverage toggle slider
                $("#log-coverage").attr("disabled", false);

                that._drawSelectedBarcode(flowcellId, txt);
            })

        }).catch(error => {
            console.error(error)
        })

    }

    /**
     * Update the master coverage chart and detail coverage chart with data from the selected barcode, triggered on change in the select box.
     * @param flowcellId {number} - the primary key of the flowcell
     * @param barcodeChosen {string} - the barcode that we are drawing the data for
     */
    _drawSelectedBarcode(flowcellId, barcodeChosen) {
        /*
        Reload a new barcode data
        */
        let that = this;
        let min, max, coverageDetail;
        let logCoverage = $("#log-coverage")[0].checked ? 1 : 0;
        let yAxisTitle = logCoverage ? "Log 10 coverage" : "Coverage";
        that._articCoverageScatterMaster.showLoading("Fetching data from server.")
        // enable zoom on the master chart
        if (!this._articCoverageScatterMaster.xAxis[0].zoomEnabled === true) {
            this._articCoverageScatterMaster.xAxis[0].zoomEnabled = true
        }
        this._axiosInstance.get('/api/v1/artic/visualisation/master', {
            params: {
                flowcellId,
                barcodeChosen,
                logCoverage
            }
        }).then(response => {
            // update the extremes on the y axis
            console.log(response.data.coverage.ymax)
            that._articCoverageScatterMaster.yAxis[0].setExtremes(0, response.data.coverage.ymax);
            that._articCoverageScatterMaster.yAxis[0].setTitle({text: yAxisTitle})
            that._updateExistingChart(that._articCoverageScatterMaster, response.data.coverage.data, 0);
            that._articCoverageScatterMaster.hideLoading()

        }).catch(error => {
            console.error(error)
        });

        coverageDetail = $("#coverageArticDetail").highcharts();
        min = coverageDetail.xAxis[0].min;
        max = coverageDetail.xAxis[0].max;
        // If min is not undefined, the charts are already displaying data, so update them, otherwise we're good
        if (min !== undefined) {
            console.log("Fetching coverage scatter dteail data");
            that._articCoverageScatterDetail.showLoading('Fetching data from the server')
            this._axiosInstance.get('/api/v1/artic/visualisation/detail', {
                params: {
                    min,
                    max,
                    barcodeChosen: that._barcodeChosen,
                    flowcellId,
                    logCoverage
                }
            }).then(response => {
                that._articCoverageScatterDetail.yAxis[0].setTitle({text: yAxisTitle})
                that._updateAxesAndData(that._articCoverageScatterDetail, min, max, response.data.coverage, true);
            }).catch(error => {
                console.error(error)
            });
        }
        this._setBarcodeMetaDataHTMLTable(flowcellId, barcodeChosen)
        this._renderPngs(flowcellId, barcodeChosen)
    }

    /**
     * @function Top level call that creates the charts on the page when the controller is initialised on page load.
     * Is called in interval to update charts.
     * @param that {ArticController} The scope of the class, allowing access to properties.
     */
    _drawArticCharts(that) {
        let min, max;
        let coverageDetail;
        let flowcellId = that._flowcellId;
        let startData;
        // If we don't have this tab, don't do all these calls.
        // if (!flowcellController.flowcellTabController.checkTabIsPresent("artic")){
        //     return
        // }
        if (!Highcharts.Series.prototype.renderCanvas) {
            throw "Module not loaded";
        }
        console.time("scatter");
        // If we are drawing it for the first time by calling the create master chart which has a callback that draws the other charts
        if (that._first === true) {
            startData = {"coverage": {"ymax": 1, "data": [0, 1]}};
            that._createCoverageMaster(startData, flowcellId);
            // Create selection drop down
            that._createOrUpdateSelectionBox(flowcellId);
            that._addListenerToLogSlider(flowcellId);
            that._drawColumnCharts()
            that._createArticPieChart(flowcellId)
            that._drawOrUpdateSummaryTable(that._coverageSummaryTable, flowcellId)
            that._first = false;
        } else {
            // update the column charts
            that._updateColumnCharts(flowcellId)
            that._updateArticPieChart(flowcellId)
            // update the table
            that._drawOrUpdateSummaryTable(that._coverageSummaryTable, flowcellId)
            // update the selection box
            that._createOrUpdateSelectionBox(flowcellId)
            //update the detail master charts
            coverageDetail = $("#coverageArticDetail").highcharts();
            // Get the minimum and maximum values of the xAxis
            min = coverageDetail.xAxis[0].min;
            // if min is undefined we don't have a detail chart so don't update
            if (this._barcodeChosen) {
                that._setBarcodeMetaDataHTMLTable(flowcellId, that._barcodeChosen)
                that._drawSelectedBarcode(flowcellId, that._barcodeChosen)
                that._renderPngs(flowcellId, that._barcodeChosen)
            }
        }
        console.timeEnd("scatter");
    }

    /**
     * @function Request the data for the column charts on the Artic Tab page and draw them.
     */
    _drawColumnCharts() {
        let that;
        let flowcellId = this._flowcellId;
        if (!this._perBarcodeCoverageChart) {
            this._perBarcodeCoverageChart = makeColumnChart(
                "artic-per-barcode-coverage",
                "PASS READ COUNTS PER BARCODE",
                "READ COUNT"
            );
        }
        if (!this._perBarcodeAverageReadLengthChart) {
            this._perBarcodeAverageReadLengthChart = makeColumnChart(
                "artic-per-barcode-average-length",
                "MEAN PASS READ LENGTH PER BARCODE",
                "READ LENGTH (BASES)"
            );
        }
        this._updateColumnCharts(flowcellId)
    }

    /**
     * Update the column charts at the top of the page.
     * @param flowcellId {number} The number of the flowcell charts
     * @private
     */
    _updateColumnCharts(flowcellId) {

        // Callback scope proof access to class methods and variables.
        this._axiosInstance.get("/api/v1/artic/visualisation/column-charts", {
            "params": {flowcellId}
        }).then(response => {
            let data, chartSeriesCoverage, chartSeriesReadLength, categories;
            let options;
            let indexInArray
            data = response.data;
            chartSeriesCoverage = this._perBarcodeCoverageChart.series;
            chartSeriesReadLength = this._perBarcodeAverageReadLengthChart.series;
            options = {
                xAxis: {categories: data["barcodes"]}
            }
            this._perBarcodeCoverageChart.update(options)
            this._perBarcodeAverageReadLengthChart.update(options)
            // we only need to check the one chart to know if both have data
            indexInArray = chartSeriesCoverage.findIndex(obj => obj.name === "Read Counts");
            // See if we already have a series by this name, if we do update existing series // show loading on the chart.
            this._perBarcodeCoverageChart.showLoading(`<div class="spinner-border text-success" role="status">
                        <span class = "sr-only"> Loading...</span></div>`);
            this._perBarcodeAverageReadLengthChart.showLoading(`<div class="spinner-border text-success" role="status">
                        <span class = "sr-only"> Loading...</span></div>`);
            if (indexInArray >= 0) {
                // if the data is not new, i.e it is identical don't draw
                if (checkHighChartsDataIsNew(data["read_counts"], chartSeriesCoverage[indexInArray].options.data)) {
                    console.log("data is identical, don't draw")
                } else {
                    chartSeriesCoverage[indexInArray].setData(data["read_counts"]);
                    chartSeriesReadLength[indexInArray].setData(data["average_read_length"])
                }
            } else { // Add the series
                this._perBarcodeCoverageChart.addSeries({"name": "Read Counts", "data": data["read_counts"]})
                this._perBarcodeAverageReadLengthChart.addSeries({
                    "name": "Average read length",
                    "data": data["average_read_length"]
                })
            }
            this._perBarcodeCoverageChart.hideLoading()
            this._perBarcodeAverageReadLengthChart.hideLoading()
        }).catch(errorResponse => {
            console.error(errorResponse)
        });
    }

    /**
     * Draw the Artic Summary Datatable, showing the Summary information for each barcode.
     * @param datatableObj {obj} DataTables object kept in the class.
     * @param flowcellId {number} The primary key of the flowcell
     * @private
     */
    _drawOrUpdateSummaryTable(datatableObj, flowcellId) {
        // If the table already exists, use the DataTable APi to update in place
        if ($.fn.DataTable.isDataTable(datatableObj)) {
            // null for callback function, false for reset paging
            datatableObj.DataTable().ajax.reload(null, false)
        } else {
            // else the datatable must be initialised
            datatableObj.DataTable({
                    // callback run on each row create
                    "createdRow": (row, data, index) => {
                        if (data.has_finished) {
                            $(row).addClass(["finished-pipeline", "artic-tr-row"]);
                        } else {
                            $(row).addClass("artic-tr-row")
                        }
                        $(row).on("click", () => {
                            this._barcodeChosen = data["barcode_name"];
                            // set value on select box
                            $("#select-artic-barcode").val(data["barcode_name"])
                            // remove placeholder option "Choose barcode"
                            $("#place-holder").remove()
                            // enable the log coverage toggle slider
                            $("#log-coverage").attr("disabled", false);
                            this._drawSelectedBarcode(this._flowcellId, data["barcode_name"])
                        })

                    },
                    "ajax": {
                        "url": "/api/v1/artic/visualisation/summary-table-data",
                        "data": {
                            flowcellId
                        },
                        async: true,
                        error: (xhr, error, code) => {
                            console.error(xhr);
                            console.error(code);
                        }
                    },
                    "columnDefs": [
                        {"targets": 0, "data": "barcode_name"},
                        {"targets": 1, "data": "chromosome__line_name"},
                        {"targets": 2, "data": "reference_line_length"},
                        {"targets": 3, "data": "read_count"},
                        {"targets": 4, "data": "total_length"},
                        {"targets": 5, "data": "average_read_length"},
                        {"targets": 6, "data": "coverage"},
                        {"targets": 7, "data": "percent_200x"},
                        {"targets": 8, "data": "percent_250x"},
                        {"targets": 9, "data": "has_sufficient_coverage"},
                    ]
                }
            );

        }
    }

    /**
     * Set the barcode metadata html table for the selected barcode
     * @param flowcellId {number} Primary key of the flowcell record in the database.
     * @param selectedBarcode {string} The buser selected barcode in the drop down.
     * @private
     */
    _setBarcodeMetaDataHTMLTable(flowcellId, selectedBarcode) {
        this._axiosInstance.get("/api/v1/artic/barcode-metadata", {
            params: {
                flowcellId,
                selectedBarcode
            }
        }).then(response => {
            $("#articMetaInfo").html(response.data)
        }).catch(error => {
            console.error(error)
        })
    }

    /**
     * Manually trigger the artic command to be run on this barcode.
     * @param flowcellId {number} The primary key of the flowcell we are running this on.
     * @param barcodeNameString {string} The primary key of the barcode that we want the data from.
     * @param event {obj} Event object created on trigger.
     */
    manuallyTriggerPipeline(flowcellId, barcodeNameString, event) {
        let jobTypeId = 17
        event.preventDefault()
        console.log(barcodeNameString.id)
        this._axiosInstance.post("/api/v1/artic/manual-trigger/", {
            flowcellId,
            barcodeName: barcodeNameString.id,
            jobTypeId
        }).then(response => {
            console.log(response)
            let message = response.data
            $("#manual-trigger").remove()
            $("#artic-message").css("display", "block").html(message)
        }).catch(error => {
            console.error(error)
        })
    }

    /**
     * Submit the user choices from the modal to the server then offer the selected files as a .tar.gz
     * @param flowcellId {number} Primary key of the flowcell
     * @param selectedBarcode {string} The chosen barcode name
     * @param event {obj} Auto generatee event object
     */
    buildResults(flowcellId, selectedBarcode, event) {
        let data = {"string": $(".results-builder-form").serialize()}
        event.preventDefault()
        data["params"] = {flowcellId, selectedBarcode}
        console.log($(".results-builder-form").serialize());
        this._axiosInstance.get("/api/v1/artic/build-results",
            {
                params: data,
                Accept: 'application/gzip',
                'Content-Type': 'application/gzip',
                responseType: 'blob'
            }).then(response => {
            const url = window.URL.createObjectURL(new Blob([response.data]));
            const link = document.createElement('a');
            link.href = url;
            link.setAttribute('download', response.headers["x-file-name"]);
            document.body.appendChild(link);
            link.click();
        }).catch(error => {
            console.error(error)
        })
    }

    /**
     * Show the barplot and boxplot PNGs generated by artic on the page it they have been made
     * @param flowcellId {number} Primary key of the flowcell
     * @param selectedBarcode {string} The chosen barcode name
     * @private
     */
    _renderPngs(flowcellId, selectedBarcode) {
        if (!$(".rendered").length && (!$("#artic-message").length || $("#marked-rerun").length)) {
            this._axiosInstance.get("/api/v1/artic/fetch-png-html", {
                params: {
                    flowcellId,
                    selectedBarcode
                }
            }).then(response => {
                console.log(response)
                $("#artic-pngs").html(response.data);
            })
        }
    }


    /**
     * Add on change listener to the HTML slider for logging coverage to call the draw selected barocde function
     * @param flowcellId {number} The primary key of the flowcell record in the database.
     * @private
     */
    _addListenerToLogSlider(flowcellId) {
        $("#log-coverage").on("change", () => {
            // set logarithmic on Y axis highcharts
            this._logged = !this._logged
            this._drawSelectedBarcode(flowcellId, this._barcodeChosen)
        })
    }

    /**
     * Rerun artic command with the data we have accrued since it's completion
     * @param flowcellId {number} The primary key of the flowcell in the database.
     * @param selectedBarcode {string} The barcode to rerun the command on.
     * @param event {Object} The event object
     */
    reRunArticCommand(flowcellId, selectedBarcode, event) {
        event.preventDefault()
        this._axiosInstance.patch("/api/v1/artic/rerun-command/", {flowcellId, selectedBarcode}).then(
            response => {
                if (!response.status === 204) {
                    console.error(`Unexpected response status. Expected 204, received ${response.status}`)
                    $("#artic-message").html("Errrorrrr.")
                    setTimeout(5000, () => {
                        $("#artic-message").empty()
                    })
                    return
                }
                $("#artic-message").html("Successfully listed artic task for rerunning.")
                $("#rerun-btn").remove()
                console.log("Successfully did a thing")
            })
    }
A
    /**
     * Create the highcharts pie chart for barcode proportion of reads unclassified vs. classified
     * @param flowcellId {number} The primary key of the flowcell in the database.
     * @private
     */
    _createArticPieChart(flowcellId) {
        if (!this._articProportionPieChart) {
            this._articProportionPieChart = makePieChart("artic-per-barcode-classification", "PROPORTION OF READS BARCODE CLASSIFIED")
            console.log("making")
        }
        this._updateArticPieChart(flowcellId)
    }

    /**
     * Update the date in highcharts pie chart for barcode proportion of reads unclassified vs. classified
     * @param flowcellId {number} The primary key of the flowcell in the database.
     * @private
     */
    _updateArticPieChart(flowcellId) {
        this._axiosInstance.get("/api/v1/artic/pie-chart-data", {params: {flowcellId}}).then(
            response => {
                let seriesData = response.data
                console.log(seriesData)
                let oldChartData = [], oldChartSeries;
                if (this._articProportionPieChart.series.length) {
                    oldChartSeries = this._articProportionPieChart.series
                    oldChartData = oldChartSeries.filter(series => {
                        return series.name === seriesData.name;
                    })
                }

                if (!oldChartData.length) {
                    this._articProportionPieChart.addSeries(seriesData[0])
                } else {
                    // if the data is the same
                    if (checkHighChartsDataIsNew(oldChartData[0].options.data, seriesData[0].data.map(e => e.y))) {
                        console.log("Skipping draw")
                    } else {
                        oldChartData[0].setData(seriesData[0].data)
                    }
                }
            }
        )
    }
}
