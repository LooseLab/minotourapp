/**
 *
 */
class ArticController {
    constructor(flowcellId) {
        this._coverageScatterMaster;
        this._coverageScatterDetail;
        this._flowcellId = flowcellId;
    }

    /**
     * Draw the master coverage chart for the selected barcode
     * @param flowcellId number - the primary key of the flowcell
     * @param barcodeChosen str - the barcode that we are drawing the data for
     */
    drawSelectedBarcode(flowcellId, barcodeChosen) {
        /*
        Reload a new chromosomes data
        */
        const axiosInstance = axios.create({
            headers: {'X-CSRFToken': csrftoken}
        });
        let self = this;
        console.log(self);
        axiosInstance.get('/api/v1/artic/coverage/master', {params: {flowcellId, barcodeChosen}}, function (data) {
            self._coverageScatterMaster.yAxis[0].setExtremes(0, data.coverage.ymax);
            updateExistingChart(self._coverageScatterMaster, data.coverage.data, 0);

        });
    }

    /**
     *
     * @param chart - Chart object that is being updated
     * @param newSeries - The new data to add into the chart
     * @param index - The index of the series that we are updating. Use 0.
     */
    updateExistingChart(chart, newSeries, index) {
        // Set the data on an existing series on a chart, either master or detail
        chart.series[index].setData(newSeries);

    }

//////////////////////////////////////////////////////////////////////////
///////// Draw the single master chart under the coverage detail /////////
//////////////////////////////////////////////////////////////////////////
    /**
     * The small coverage master chart that shows the region selected
     * @param data
     * @param flowcellId
     */
    createCoverageMaster(data, flowcellId) {
        const self = this;
        this._coverageScatterMaster = Highcharts.chart("coverageScatterMaster", {

            chart: {
                zoomType: "x",
                height: "150px",
                marginLeft: 50,
                marginRight: 20,
                events: {
                    selection: afterSelection.bind(this)
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
        }, (self) => {
            ////////////////////////////////////////////////////////////
            ///////// Call the detail chart drawing functions //////////
            ////////////////////////////////////////////////////////////
            coverageScatterDetail = self.createDetailChart("coverageScatterDetail", data.coverage.ymax, "Coverage");
        });

    }

    /**
     *
     * @param targetDivId str Div ID chart is drawn to
     * @param ymax int The maximum value of the Y axis
     * @param title str The chart title
     * @param chartType str Chart title
     * @param stepped bool Stepped line
     * @param legend bool Display chart legend
     * @returns {*}
     */
    createDetailChart(targetDivId, ymax, title, chartType = "scatter", stepped = false, legend = false) {
        // The data we are going to use to display, starts a display with an empty chart
        let detailData = [];
        // The variable assigned to the chart we will return
        let chartyChart;

        chartyChart = Highcharts.chart(targetDivId, {

            chart: {
                zoomType: "x",
                height: "300px",
                reflow: false,
                events: {
                    selection: afterSelection.bind(this)
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
                    return 'The value for base <b>' + this.x +
                        '</b> is <b>' + this.y + '</b>';
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
        let flowcellId = this._flowcellId;
        // Prevent the default behaviour
        event.preventDefault();
        // The minimum and maximum value on the X axis
        var min = Math.trunc(event.xAxis[0].min);
        var max = Math.trunc(event.xAxis[0].max);
        // If the selected box is too large tell the selector to learn to be less greedy
        if (max - min > showWidth) {
            label = coverageScatterDetail.renderer.label('Please select a smaller region', 100, 120)
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
        coverageScatterMaster.xAxis[0].removePlotBand('mask-before');
        // Add the new one in the correct location
        coverageScatterMaster.xAxis[0].addPlotBand({
            id: 'mask-before',
            from: min,
            to: max,
            color: 'rgba(0, 0, 0, 0.2)'
        });

        // Show the loading screen on the detail charts
        let charts = [coverageScatterDetail];

        charts.forEach(function (chart) {
            chart.showLoading('Fetching data from the server');
        });

        // Get the detail chart data for the zoom selection from the server
        $.getJSON('/api/v1/artic/coverage/detail', {
            min,
            max,
            chromosome_chosen,
            flowcellId
        }, function (newSeriesDict, statusText, xhr) {

            // if successful request, set the detail charts to show the newly retrieve data
            if (xhr.status === 200) {

                update_axes_and_data(coverageScatterDetail, min, max, newSeriesDict.coverage, true);

            } else if (xhr.status === 404) {

            } else {
                throw "Request error";
            }

        });


    }

    // draw the charts, this is the function that is called in request data every 60 seconds
    drawArticChart() {
        let flowcellId = document.querySelector("#flowcell-id").value;
        create_eb_selection_box(flowcellId);

        if (!Highcharts.Series.prototype.renderCanvas) {
            throw "Module not loaded";
        }

        // Select the coverage detail chart to see if it already exists
        let coverageDetail = $("#coverageScatterDetail").highcharts();

        let min, max;
        console.time("scatter");

        let startData = {"coverage": {"ymax": 1, "data": [0, 1]}};
        // If the coverage detail chart doesn't exist, draw it for the first time by calling the create master chart which has a callback that draws the other charts
        if (coverageDetail === undefined) {
            createCoverageMaster(startData, flowcellId);
        } else {
            // If the coverage detail chart does exist, update the existing charts
            // Get the minimum and maximum values of the xAxis
            min = coverageDetail.min;
            max = coverageDetail.max;
            // If min is not undefined, the charts are already displaying data, so update thew, otherwise we're good
            if (min !== undefined) {
                $.getJSON('/api/v1/artic/coverage/detail', {
                    min,
                    max,
                    chromosome_chosen,
                    flowcellId
                }, function (newSeriesDict, statusText, xhr) {
                    updateExistingChart(coverageScatterMaster, data.coverage.data, 0);
                    updateExistingChart(coverageDetail, newSeriesDict.coverage.data, 0);
                });
            }


        }
        console.timeEnd("scatter");

    }

}
