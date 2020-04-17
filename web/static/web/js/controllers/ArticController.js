/**
 *
 */
class ArticController {
    constructor(divName) {
        this._coverageScatterMaster;
        this._coverageScatterDetail;
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

        axiosInstance.get('/api/v1/artic/coverage/master', {params: {flowcellId, barcodeChosen}}, function (data) {
            coverageScatterMaster.yAxis[0].setExtremes(0, data.coverage.ymax);
            updateExistingChart(coverageScatterMaster, data.coverage.data, 0);

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
        this._flowcellId = flowcellId;
        coverageScatterMaster = Highcharts.chart("coverageScatterMaster", {

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
        }, function (masterChart) {
            ////////////////////////////////////////////////////////////
            ///////// Call the detail chart drawing functions //////////
            ////////////////////////////////////////////////////////////
            coverageScatterDetail = createDetailChart("coverageScatterDetail", data.coverage.ymax, "Coverage");
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

}
