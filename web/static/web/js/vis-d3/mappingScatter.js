// Initialise the plots so they are global variables
let coverageScatterDetail;
let coverageScatterMaster;
let matchScatterDetail;
let mismatchScatterDetail;
let localBenefitDetail;
let rollingBenefitLineDetail;
let forwardMaskDetail;
let revMaskDetail;

// showWidth is the maximum number of bases we can select to view on the master chart
let showWidth = 250000;

function updateExistingChart(chart, newSeries, index) {
    // Set the data on an existing series on a chart, either master or detail
    chart.series[index].setData(newSeries);
}

////////////////////////////////////////////////////////////////////
///////// Draw the detail charts under the coverage detail /////////
////////////////////////////////////////////////////////////////////

function createChart(targetDivId, ymax, title, chartType = "scatter", stepped = false, legend=false) {
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
            min: 0,
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

//////////////////////////////////////////////////////////////////////////
///////// Draw the single master chart under the coverage detail /////////
//////////////////////////////////////////////////////////////////////////
function createCoverageMaster(data) {
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
        coverageScatterDetail = createChart("coverageScatterDetail", data.coverage.ymax, "Coverage");
        matchScatterDetail = createChart("matchScatterDetail", data.match.ymax, "Matches");
        mismatchScatterDetail = createChart("mismatchScatterDetail", data.mismatch.ymax, "Mismatches");
        localBenefitDetail = createChart("localBenScatterDetail", data.localBen.ymax, "Local benefit");
        rollingBenefitLineDetail = createChart("rollingBenefitLineDetail", data.forwardRoll.ymax, "Rolling benefit", "line", false, true);
        forwardMaskDetail = createChart("maskStepLineFwdDetail", data.fwdMask.ymax, "Forward mask", "line", true);
        revMaskDetail = createChart("maskStepLineRevDetail", data.revMask.ymax, "Reverse mask", "line", true);
    });
}

function afterSelection(event) {
    // Event function called after the zoom box is chosen
    let label;
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
    coverageScatterDetail.showLoading('Fetching data from the server');
    matchScatterDetail.showLoading('Fetching data from the server');
    mismatchScatterDetail.showLoading('Fetching data from the server');
    localBenefitDetail.showLoading('Fetching data from the server');
    rollingBenefitLineDetail.showLoading('Fetching data from the server');
    forwardMaskDetail.showLoading('Fetching data from the server');
    revMaskDetail.showLoading('Fetching data from the server');

    // Get the detail chart data for the zoom selection from the server
    $.getJSON('/api/v1/readuntil/benefitdata/detail', {min, max}, function (newSeriesDict, statusText, xhr) {

        // if successful request, set the detail charts to show the newly retrieve data
        if (xhr.status === 200) {
            coverageScatterDetail.series[0].setData(newSeriesDict.coverage.data);
            coverageScatterDetail.min = min;
            coverageScatterDetail.max = max;
            coverageScatterDetail.hideLoading();

            matchScatterDetail.series[0].setData(newSeriesDict.match.data);
            matchScatterDetail.min = min;
            matchScatterDetail.max = max;
            matchScatterDetail.hideLoading();

            mismatchScatterDetail.series[0].setData(newSeriesDict.mismatch.data);
            mismatchScatterDetail.min = min;
            mismatchScatterDetail.max = max;
            mismatchScatterDetail.hideLoading();

            localBenefitDetail.series[0].setData(newSeriesDict.localBenefit.data);
            localBenefitDetail.min = min;
            localBenefitDetail.max = max;
            localBenefitDetail.hideLoading();

            rollingBenefitLineDetail.series[0].setData(newSeriesDict.rollingBenefitFwd.data);
            rollingBenefitLineDetail.series[0].name = "Forward rolling benefit";
            rollingBenefitLineDetail.min = min;
            rollingBenefitLineDetail.max = max;
            rollingBenefitLineDetail.hideLoading();

            // check if we have a reverse line already or not on the rolling benefit detail chart
            if (rollingBenefitLineDetail.series[1] === undefined){
                // if there is no line, add a new series to the forward reverse rolling benefit line for the reverse line data
                rollingBenefitLineDetail.addSeries({name: "Reverse rolling benefit",
                data:newSeriesDict.rollingBenefitRev.data, "color": "gold"});
            } else {
                // if there is a line, update the already existing series on the forward line chart
                rollingBenefitLineDetail.series[1].setData(newSeriesDict.rollingBenefitRev.data);
            }


            forwardMaskDetail.series[0].setData(newSeriesDict.forwardMask.data);
            forwardMaskDetail.min = min;
            forwardMaskDetail.max = max;
            forwardMaskDetail.hideLoading();

            revMaskDetail.series[0].setData(newSeriesDict.reverseMask.data);
            revMaskDetail.min = min;
            revMaskDetail.max = max;
            revMaskDetail.hideLoading();

        } else {
            throw "Request error";
        }
    });


}
// draw the charts, this is the function that is called in request data every 60 seconds
function drawReadUntilCharts() {
    $.getJSON('/api/v1/readuntil/benefitdata/master', function (data) {
        if (!Highcharts.Series.prototype.renderCanvas) {
            throw "Module not loaded";
        }

        // Select the coverage detail chart to see if it already exists
        let coverageDetail = $("#coverageScatterDetail").highcharts();
        let min, max;
        // Wake up, Neo....
        // The matrix has you
        // Follow the white rabbit
        console.time("scatter");
        // If the coverage detail chart doesn't exist, draw it for the first time by calling the create master chart which has a callback that draws the other charts
        if (coverageDetail === undefined) {
            createCoverageMaster(data);
        } else {
            // If the coverage detail chart does exist, update the existing charts
            // Get the minimum and maximum values of the xAxis
            min = coverageDetail.min;
            max = coverageDetail.max;
            // If min is not undefined, the charts are already displaying data, so update the, otherwise we're good
            if (min !== undefined){
                $.getJSON('/api/v1/readuntil/benefitdata/detail', {min, max}, function (newSeriesDict, statusText, xhr) {
                    updateExistingChart(coverageScatterMaster, data.coverage.data, 0);
                    updateExistingChart(coverageDetail, newSeriesDict.coverage.data, 0);
                    updateExistingChart(matchScatterDetail, newSeriesDict.match.data, 0);
                    updateExistingChart(mismatchScatterDetail, newSeriesDict.mismatch.data, 0);
                    updateExistingChart(localBenefitDetail, newSeriesDict.localBenefit.data, 0);
                    updateExistingChart(rollingBenefitLineDetail, newSeriesDict.rollingBenefitFwd.data, 0);
                    updateExistingChart(rollingBenefitLineDetail, newSeriesDict.rollingBenefitRev.data, 1);
                    updateExistingChart(forwardMaskDetail, newSeriesDict.forwardMask.data, 0);
                    updateExistingChart(revMaskDetail, newSeriesDict.reverseMask.data, 0);
            });
            }


        }
        console.timeEnd("scatter");

    });
}