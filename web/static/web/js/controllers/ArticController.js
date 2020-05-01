/**
 *
 */
class ArticController {
    constructor(flowcellId) {
        console.log("Initialising Artic controoller");
        this._coverageScatterMaster;
        this._coverageScatterDetail;
        this._flowcellId = flowcellId;
        // this.createCoverageMaster();
        this.drawArticChart(flowcellId);
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
        let that = this;
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
        }, () => {
            ////////////////////////////////////////////////////////////
            ///////// Call the detail chart drawing functions //////////
            ////////////////////////////////////////////////////////////
            coverageScatterDetail = that.createDetailChart("coverageScatterDetail", data.coverage.ymax, "Coverage");
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

    createEbSelectionBox(flowcellId) {
        let that = this
        // For all the sel objects that finding the element by the sel class returned
        $('.sel.artic').each(function () {
            // get the child of the div, which is the select option in this case, do not display it
            $(this).children('select').css('display', 'none');
            // Get the sel div element as a var
            var $current = $(this);
            // get all the options children
            $.get("/api/v1/artic/visualisation/barcodes", {flowcellId}, function (barcodes, statusText, xhr) {
                // if it's the first options clause
                if (xhr.status != 200) console.error("Error");
                // If we're on our first drawing of the select spans
                let first_iteration = !document.getElementById("artic-placeholder");
                // The barcodes we have that are already present
                let alreadyPresentBarcode = [];
                // If it is not our first drawing of it
                if (!first_iteration) {
                    // The already present barcodes html elements
                    alreadyPresentBarcode = [...$(".sel__box")[0].children];
                    // for each of the elements
                    alreadyPresentBarcode.forEach(function (presBarcode) {
                        // if the chromsome in already present chromsome is in the list of barcodes fetched from the server
                        if (barcodes.includes(presBarcode.textContent)) {
                            // get the index,
                            let index = barcodes.indexOf(presBarcode.textContent);
                            // remove it from the freshly fetched barcodes
                            barcodes.splice(index, 1);
                        }
                    });
                    // If there are no new barcodes
                    if (!barcodes.length) {
                        return;
                    }
                    // if it is our first iteration
                } else {
                    // unshift choose barcodes
                    barcodes.unshift("Choose Barcode");
                }


                barcodes.forEach(function (barcode, index) {
                    console.log(barcode)
                    // if this is the first time so we're creating the arrays
                    if (index == 0 && first_iteration) {
                        // the div to the front of the sel div, make it sel__box class instead of sel
                        $current.prepend($('<div>', {
                            class: $current.attr('class').replace(/sel/g, 'sel__box ')

                        }));
                        // Set the placeholder text to the option text - This would be the first barcode
                        var placeholder = barcode;
                        // prepend a span in front of the box div, with the placeholder text
                        $current.prepend($('<span>', {
                            class: $current.attr('class').replace(/sel/g, 'sel__placeholder '),
                            text: placeholder,
                            id: "artic-placeholder",
                            'data-placeholder': placeholder
                        }));

                        return;
                    }

                    // Then for the rest of the options add them in
                    $current.children('div').append($('<span>', {
                        class: $current.attr('class').replace(/sel/g, 'sel__box__options '),
                        text: barcode
                    }));
                });
                // Toggling the `.active` state on the `.sel`.
                $('.sel').click(function () {
                    $(this).toggleClass('active');
                });

                // Toggling the `.selected` state on the options.
                $('.sel__box__options').click(function () {
                    console.log("The things....")
                    var txt = $(this).text();
                    chromosome_chosen = txt;
                    var index = $(this).index();


                    $(this).siblings('.sel__box__options').removeClass('selected');
                    $(this).addClass('selected');

                    var $currentSel = $(this).closest('.sel');
                    $currentSel.children('.sel__placeholder').text(txt);
                    $currentSel.children('select').prop('selectedIndex', index + 1);

                    that.drawSelectedBarcode(flowcellId, txt);


                });
                $(".sel").css("display", "");
            });
        });
    }

    /**
     * Draw the master coverage chart for the selected barcode
     * @param flowcellId number - the primary key of the flowcell
     * @param barcodeChosen str - the barcode that we are drawing the data for
     */
    drawSelectedBarcode(flowcellId, barcodeChosen){
        /*
        Reload a new chromosomes data
        */
        let that = this;
        console.log(that._coverageScatterMaster);
        $.getJSON('/api/v1/artic/visualisation/master', {flowcellId, barcodeChosen}, function (data) {
            that._coverageScatterMaster.yAxis[0].setExtremes(0, data.coverage.ymax);
            that.updateExistingChart(that.coverageScatterMaster, data.coverage.data, 0);

        });
    }

    // draw the charts, this is the function that is called in request data every 60 seconds
    drawArticChart() {
        let that = this;
        let flowcellId = document.querySelector("#flowcell-id").value;
        that.createEbSelectionBox(flowcellId);

        if (!Highcharts.Series.prototype.renderCanvas) {
            throw "Module not loaded";
        }

        // Select the coverage detail chart to see if it already exists
        let coverageDetail = $("#coverageArticDetail");

        let min, max;
        console.time("scatter");

        let startData = {"coverage": {"ymax": 1, "data": [0, 1]}};
        // If the coverage detail chart doesn't exist, draw it for the first time by calling the create master chart which has a callback that draws the other charts
        if (coverageDetail.length === 0) {
            that.createCoverageMaster(startData, flowcellId);
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
                    that.updateExistingChart(that.coverageScatterMaster, data.coverage.data, 0);
                    that.updateExistingChart(coverageDetail, newSeriesDict.coverage.data, 0);
                });
            }


        }
        console.timeEnd("scatter");

    }

}
