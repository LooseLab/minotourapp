/**
 *
 */
class ArticController {
    constructor(flowcellId) {
        this._articCoverageScatterMaster = "hello";
        this._articCoverageScatterDetail = "world";
        this._flowcellId = flowcellId;
        this._first = true;
        this._showWidth = 250000;
        this._barcodeChosen = "banter"
        this._axiosInstance = axios.create({
            headers: {"X-CSRFToken": getCookie('csrftoken')}
        });
        this._barcodeSelect = $("#select-artic-barcode");
        this._perBarcodeAverageReadLengthChart = null;
        this._perBarcodeCoverageChart = null;
        this._coverageSummaryTable = $("#artic-coverage-summary-table");
        // this.createCoverageMaster();
        this.drawArticChart(flowcellId);
    }


    /**
     *
     * @param {obj} chart - Chart object that is being updated
     * @param {obj} newSeries - The new data to add into the chart
     * @param {number} index - The index of the series that we are updating. Use 0.
     */
    updateExistingChart(chart, newSeries, index) {
        // Set the data on an existing series on a chart, either master or detail
        chart.series[index].setData(newSeries);

    }


    /**
     * The small coverage master chart that shows the region selected
     * @param data
     * @param flowcellId
     */
    createCoverageMaster(data, flowcellId) {
        const self = this;
        let that = this;
        this._articCoverageScatterMaster = Highcharts.chart("coverageArticMaster", {

            chart: {
                zoomType: "x",
                height: "150px",
                marginLeft: 50,
                marginRight: 20,
                events: {
                    selection: that.afterSelection.bind(this)
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
                maxZoom: 10000, // fourteen days
                plotBands: [{
                    id: 'mask-before',
                    from: data.coverage.data[0][0],
                    to: data.coverage.data[data.coverage.data.length - 1][0],
                    color: 'rgba(0, 0, 0, 0.2)'
                }]
            },

            yAxis: {
                // Renders faster when we don"t have to compute min and max
                min: 0,
                max: data.coverage.ymax,
                title: {
                    text: null
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
            // Call back function, called after the master chart is drawn, draws all the detail charts
        }, () => {

            that._articCoverageScatterDetail = that.createDetailChart("coverageArticDetail", data.coverage.ymax, "Coverage", "line"
                , false);
        });

    }

    /**
     *
     * @param {string} targetDivId Div ID chart is drawn to
     * @param {number} ymax The maximum value of the Y axis
     * @param {string} title The chart title
     * @param {string} chartType Chart title
     * @param {boolean} stepped Stepped line
     * @param {boolean} legend Display chart legend
     * @returns {obj} chartyChart Highcharts chart.
     */
    createDetailChart(targetDivId, ymax, title, chartType = "scatter", stepped = false, legend = false) {
        // The data we are going to use to display, starts a display with an empty chart
        let detailData = [];
        // The variable assigned to the chart we will return
        let chartyChart;

        let that = this;

        chartyChart = Highcharts.chart(targetDivId, {

            chart: {
                zoomType: "x",
                height: "300px",
                reflow: false,
                events: {
                    selection: that.afterSelection.bind(this)
                }
            },

            boost: {
                useGPUTranslations: true,
                usePreAllocated: true
            },

            xAxis: {
                type: chartType,
                gridLineWidth: 1,
                text: {
                    text: "Base number"
                }
            },

            yAxis: {
                // Renders faster when we don"t have to compute min and max
                softMin: 0,
                max: ymax,
                title: {
                    text: null
                }
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
            exporting: {"enabled": false},

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

    afterSelection(event) {
        // Event function called after the zoom box is chosen
        let label;
        let that = this;
        let flowcellId = this._flowcellId;
        let charts = [this._articCoverageScatterDetail];
        // The minimum and maximum value on the X axis
        let min = Math.trunc(event.xAxis[0].min);
        let max = Math.trunc(event.xAxis[0].max);
        // Prevent the default behaviour
        event.preventDefault();

        // If the selected box is too large tell the selector to learn to be less greedy
        if (max - min > that._showWidth) {
            label = that._articCoverageScatterDetail.renderer.label('Please select a smaller region', 100, 120)
                .attr({
                    fill: Highcharts.getOptions().colors[0],
                    padding: 10,
                    r: 5,
                    zIndex: 8
                })
                .css({
                    color: '#FFFFFF'
                })
                .add();

            setTimeout(function () {
                label.fadeOut();
            }, 2000);
            return;
        }

        // Change the grey plot band on the master chart to reflect the selected region
        that._articCoverageScatterMaster.xAxis[0].removePlotBand('mask-before');
        // Add the new one in the correct location
        that._articCoverageScatterMaster.xAxis[0].addPlotBand({
            id: 'mask-before',
            from: min,
            to: max,
            color: 'rgba(0, 0, 0, 0.2)'
        });

        // Show the loading screen on the detail charts

        charts.forEach(function (chart) {
            chart.showLoading('Fetching data from the server');
        });

        this._axiosInstance.get('/api/v1/artic/visualisation/detail', {
            params: {min, max, barcodeChosen: this._barcodeChosen, flowcellId}
        }).then(response => {
            // if successful request, set the detail charts to show the newly retrieve data
            if (response.status === 200) {

                that._updateAxesAndData(that._articCoverageScatterDetail, min, max, response.data.coverage, true);

            } else if (response.status === 404) {

            } else {
                throw "Request error";
            }

        }).catch(error => {
            console.error(error)
        });


    }


    _updateAxesAndData(chart, min, max, data, ymax) {
        // Y max let's us know we need to change the yAxis
        // Set the xAxis extremes to the newly selected values
        chart.xAxis[0].setExtremes(min, max);
        // if we need to update the yaxis as well
        if (ymax) {
            // set the yAxisextremes
            chart.yAxis[0].setExtremes(0, data.ymax);
        }
        // set the data for the chart to display
        chart.series[0].setData(data.data);
        // update the axes
        chart.xAxis[0].update();

        // hide the loading screen
        chart.hideLoading();
    }

    /**
     * Draw the select box on the artic page for the coverage tracker.
     * @param flowcellId {number} Primary key of the flowcell that barcodes are bring added to the selection box for.
     * @private
     */
    _createSelectionBox(flowcellId) {
        let first_iteration, optionChildren, barcodes, alreadyPresentBarcode, that, $current, barcodeList;
        that = this;
        barcodeList = [];
        // for each select box on the page
        $current = $("#select-artic-barcode");
        optionChildren = $current.children('option');
        first_iteration = optionChildren.length === 0 ? true : false;
        that._axiosInstance.get("/api/v1/artic/visualisation/barcodes", {params: {flowcellId}}).then((response) => {
            if (response.status !== 200) console.error("Error");
            barcodes = response.data;
            // else we need to remove the barcodes that we have already added then add them.
            if (!first_iteration) {
                alreadyPresentBarcode = [...$("#select-artic-barcode")[0].children];
                alreadyPresentBarcode.forEach(function (presBarcode) {
                    // if the chromsome in already present chromsome is in the list of barcodes fetched from the server
                    if (barcodes.includes(presBarcode.textContent)) {
                        // get the index,
                        let index = barcodes.indexOf(presBarcode.textContent);
                        // remove it from the freshly fetched barcodes
                        barcodes.splice(index, 1);
                    }
                });
            } else {
                barcodes.unshift("Choose Barcode");
            }


            barcodes.forEach(barcode => {
                barcodeList.push(`<option class="sel__box__options" value="${barcode}">${barcode}</option>`)
            })
            $current.html(barcodeList.join(""))
            $current.change(function()  {
                let txt;
                txt = $(this).val();
                that._barcodeChosen = txt;

                that.drawSelectedBarcode(flowcellId, txt);
            })

        }).catch(error => {
            console.error(error)
        })

    }
        ;


// _createSelectionBox(flowcellId) {
//     let that = this
//     // For all the sel objects that finding the element by the sel class returned
//     $('.sel.artic').each(function () {
//         // Get the sel div element as a var
//         let $current = $(this);
//         // get the child of the div, which is the select option in this case, do not display it
//         $(this).children('select').css('display', 'none');
//         // get all the options children
//         that._axiosInstance.get("/api/v1/artic/visualisation/barcodes", {params: {flowcellId}}).then((response) => {
//             let barcodes = response.data;
//             // if it's the first options clause
//             if (response.status != 200) console.error("Error");
//             // If we're on our first drawing of the select spans
//             let first_iteration = !document.getElementById("artic-placeholder");
//             // The barcodes we have that are already present
//             let alreadyPresentBarcode = [];
//             // If it is not our first drawing of it
//             if (!first_iteration) {
//                 // The already present barcodes html elements
//                 alreadyPresentBarcode = [...$(".sel__box")[0].children];
//                 // for each of the elements
//                 alreadyPresentBarcode.forEach(function (presBarcode) {
//                     // if the chromsome in already present chromsome is in the list of barcodes fetched from the server
//                     if (barcodes.includes(presBarcode.textContent)) {
//                         // get the index,
//                         let index = barcodes.indexOf(presBarcode.textContent);
//                         // remove it from the freshly fetched barcodes
//                         barcodes.splice(index, 1);
//                     }
//                 });
//                 // If there are no new barcodes
//                 if (!barcodes.length) {
//                     return;
//                 }
//                 // if it is our first iteration
//             } else {
//                 // unshift choose barcodes
//                 barcodes.unshift("Choose Barcode");
//             }
//
//
//             barcodes.forEach(function (barcode, index) {
//                 // if this is the first time so we're creating the arrays
//                 if (index == 0 && first_iteration) {
//                     // the div to the front of the sel div, make it sel__box class instead of sel
//                     $current.prepend($('<div>', {
//                         class: $current.attr('class').replace(/sel/g, 'sel__box ')
//
//                     }));
//                     // Set the placeholder text to the option text - This would be the first barcode
//                     var placeholder = barcode;
//                     // prepend a span in front of the box div, with the placeholder text
//                     $current.prepend($('<span>', {
//                         class: $current.attr('class').replace(/sel/g, 'sel__placeholder '),
//                         text: placeholder,
//                         id: "artic-placeholder",
//                         'data-placeholder': placeholder
//                     }));
//
//                     return;
//                 }
//
//                 // Then for the rest of the options add them in
//                 $current.children('div').append($('<span>', {
//                     class: $current.attr('class').replace(/sel/g, 'sel__box__options '),
//                     text: barcode
//                 }));
//             });
//             // Toggling the `.active` state on the `.sel`.
//             $('.sel').click(function () {
//                 $(this).toggleClass('active');
//             });
//
//             // Toggling the `.selected` state on the options.
//             $('.sel__box__options').click(function () {
//                 console.log("The things....")
//                 var txt = $(this).text();
//                 that._barcodeChosen = txt;
//                 var index = $(this).index();
//
//
//                 $(this).siblings('.sel__box__options').removeClass('selected');
//                 $(this).addClass('selected');
//
//                 var $currentSel = $(this).closest('.sel');
//                 $currentSel.children('.sel__placeholder').text(txt);
//                 $currentSel.children('select').prop('selectedIndex', index + 1);
//
//                 that.drawSelectedBarcode(flowcellId, txt);
//
//
//             });
//             $(".sel").css("display", "");
//         }).catch(error => {
//             console.error(error)
//         });
//     });
// }

    /**
     * Draw the master coverage chart for the selected barcode
     * @param flowcellId number - the primary key of the flowcell
     * @param barcodeChosen str - the barcode that we are drawing the data for
     */
    drawSelectedBarcode(flowcellId, barcodeChosen) {
        /*
        Reload a new chromosomes data
        */
        let that = this;
        let min, max, coverageDetail;
        this._axiosInstance.get('/api/v1/artic/visualisation/master', {
            params: {
                flowcellId,
                barcodeChosen
            }
        }).then(response => {
            that._articCoverageScatterMaster.yAxis[0].setExtremes(0, response.data.coverage.ymax);
            that.updateExistingChart(that._articCoverageScatterMaster, response.data.coverage.data, 0);

        }).catch(error => {
            console.error(error)
        });

        coverageDetail = $("#coverageArticDetail").highcharts();
        min = coverageDetail.xAxis[0].min;
        max = coverageDetail.xAxis[0].max;
        // If min is not undefined, the charts are already displaying data, so update them, otherwise we're good
        if (min !== undefined) {
            that._articCoverageScatterDetail.showLoading('Fetching data from the server')

            this._axiosInstance.get('/api/v1/artic/visualisation/detail', {
                params: {
                    min,
                    max,
                    barcodeChosen: that._barcodeChosen,
                    flowcellId
                }
            }).then(response => {
                that._updateAxesAndData(that._articCoverageScatterDetail, min, max, response.data.coverage, true);
                // that.updateExistingChart(coverageDetail, response.coverage.data, 0);
                // that.updateExistingChart(that._articCoverageScatterMaster, response.coverage.data, 0);
            }).catch(error => {
                console.error(error)
            });
        }
    }

    /**
     * @function Top level call that creates the chart om the page when the page is laoded
     */
    drawArticChart() {
        let min, max;
        let that = this;
        let coverageDetail;

        let flowcellId = this._flowcellId;
        let startData = {"coverage": {"ymax": 1, "data": [0, 1]}};


        that._createSelectionBox(flowcellId);

        if (!Highcharts.Series.prototype.renderCanvas) {
            throw "Module not loaded";
        }

        // Select the coverage detail chart to see if it already exists

        console.time("scatter");

        // If the coverage detail chart doesn't exist, draw it for the first time by calling the create master chart which has a callback that draws the other charts
        if (this._first === true) {
            that.createCoverageMaster(startData, flowcellId);
            this._first = false;
            that.drawColumnCharts()
            that._drawSummaryTable()
        } else {

            coverageDetail = $("#coverageArticDetail").highchart();
            // If the coverage detail chart does exist, update the existing charts
            // Get the minimum and maximum values of the xAxis
            min = coverageDetail.xAxis[0].min;
            max = coverageDetail.xAxis[0].max;
            // If min is not undefined, the charts are already displaying data, so update them, otherwise we're good
            if (min !== undefined) {
                this._axiosInstance.get('/api/v1/artic/visualisation/detail', {
                    params: {
                        min,
                        max,
                        barcodeChosen: that._barcodeChosen,
                        flowcellId
                    }
                }).then(response => {
                    that.updateExistingChart(coverageDetail, response.coverage.data, 0);
                    that.updateExistingChart(that._articCoverageScatterMaster, response.coverage.data, 0);
                }).catch(error => {
                    console.error(error)
                });
                this.drawSelectedBarcode(flowcellId, this._barcodeChosen)
            }


        }
        console.timeEnd("scatter");

    }

    /**
     * @function Request the data for the column charts on the Artic Tab page and draw them.
     */
    drawColumnCharts() {
        let that;
        let flowcellId = this._flowcellId;
        if (!this._perBarcodeCoverageChart) {

            this._perBarcodeCoverageChart = makeColumnChart(
                "artic-per-barcode-coverage",
                "Coverage Per Barcode",
                "Coverage X"
            );

        }

        if (!this._perBarcodeAverageReadLengthChart) {

            this._perBarcodeAverageReadLengthChart = makeColumnChart(
                "artic-per-barcode-average-length",
                "Mean Read Length By Barcode",
                "Mean Read Length By Barcode"
            );

        }
        // show loading on the page.
        this._perBarcodeCoverageChart.showLoading(`<div class="spinner-border text-success" role="status">
    <span class = "sr-only"> Loading...</span></div>`);
        this._perBarcodeAverageReadLengthChart.showLoading(`<div class="spinner-border text-success" role="status">
    <span class = "sr-only"> Loading...</span></div>`);
        // Callback scope proof access to class methods and variables.
        that = this;

        this._axiosInstance.get("/api/v1/artic/visualisation/column-charts", {
            "params": {flowcellId}
        }).then(response => {
            let data, chartSeriesCoverage, chartSeriesReadLength;


            data = response.data;

            chartSeriesCoverage = this._perBarcodeCoverageChart.series;

            chartSeriesReadLength = this._perBarcodeAverageReadLengthChart.series;

            Object.entries(data).forEach(([key, value], index) => {
                let indexInArray = chartSeriesCoverage.findIndex(obj => obj.name === key);
                // See if we already have a series by this name, if we do update exisiting series
                if (indexInArray >= 0) {
                    chartSeriesCoverage[indexInArray].setData(value["coverage"]);
                    chartSeriesCoverage[indexInArray].update({
                        name: key
                    });
                    chartSeriesReadLength[indexInArray].setData(value["average_read_length"])
                    chartSeriesReadLength[indexInArray].update({
                        name: key
                    });

                } else { // Add the series
                    this._perBarcodeCoverageChart.addSeries({"name": key, "data": [value["coverage"]]})
                    this._perBarcodeAverageReadLengthChart.addSeries({
                        "name": key,
                        "data": [value["average_read_length"]]
                    })
                }
            })
            this._perBarcodeCoverageChart.hideLoading()
            this._perBarcodeAverageReadLengthChart.hideLoading()
        }).catch(errorResponse => {
            console.error(errorResponse.message)
        })
    }


    /**
     *
     * @param data
     * @private
     */
    _updateColumnCharts(data) {

    }

    _drawSummaryTable() {
        let that = this;
        // If the table already exists, use the DataTable APi to update in place
        if ($.fn.DataTable.isDataTable(that._coverageSummaryTable)) {
            that._coverageSummaryTable.DataTable().clear();
            that._coverageSummaryTable.DataTable().draw();
        } else {
            // else the databale must be initialised
            that._coverageSummaryTable.DataTable({
                    "ajax": {
                        "url": "/api/v1/artic/visualisation/summary-table-data",
                        "data": {
                            "flowcellId": that._flowcellId
                        },
                        async: true,
                        error: (xhr, error, code) => {
                            console.error(xhr);
                            console.error(code);
                        }
                    },
                    "columns": [
                        {"data": "barcode_name"},
                        {"data": "chromosome__line_name"},
                        {"data": "reference_line_length"},
                        {"data": "read_count"},
                        {"data": "total_length"},
                        {"data": "average_read_length"},
                        {"data": "coverage"},
                    ]
                }
            );
        }
    }

}
