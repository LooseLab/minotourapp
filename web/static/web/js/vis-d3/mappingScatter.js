let coverageScatterDetail;
let coverageScatterMaster;
let matchScatterDetail;
let mismatchScatterDetail;
let localBenefitDetail;
let rollingBenefitLineDetail;
let forwardMaskDetail;
let revMaskDetail;

let showWidth = 250000;

function updateExistingChart(chart, new_series) {

    chart.series[0].setData(new_series);
}

function createChart(targetDivId, ymax) {
    let detailData = [];
    let chartyChart;
    chartyChart = Highcharts.chart(targetDivId, {

        chart: {
            zoomType: "x",
            height: "400px",
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
            type: "scatter",
            gridLineWidth: 1
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
            text: "Wake up, Neo..."
        },

        legend: {
            enabled: false
        },

        exporting: {"enabled": false},

        series: [{
            type: "scatter",
            color: "yellow",
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

function createMatchScatterDetail(masterChart, data) {
    /////////////////////////////////////////
    /// Draw the matches scatter plot ///////
    /////////////////////////////////////////
    let detailData = [];

    matchScatterDetail = Highcharts.chart("matchScatterDetail", {

        chart: {
            zoomType: "x",
            height: "400px",
            reflow: false,
            className: "coverageDetail",
            events: {
                selection: afterSelection.bind(this)
            }
        },

        boost: {
            useGPUTranslations: true,
            usePreAllocated: true
        },

        xAxis: {
            type: "scatter",
            gridLineWidth: 1
        },

        yAxis: {
            // Renders faster when we don"t have to compute min and max
            min: 0,
            max: data.match.ymax,
            title: {
                text: null
            }
        },

        title: {
            text: "The Matrix has you..."
        },

        legend: {
            enabled: false
        },

        series: [{
            type: "scatter",
            color: "green",
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
}

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
            color: "yellow",
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
    }, function (masterChart) {
        // createCoverageDetail(masterChart, data);
        // createMatchScatterDetail(masterChart, data);
        coverageScatterDetail = createChart("coverageScatterDetail", data.coverage.ymax);
        matchScatterDetail = createChart("matchScatterDetail", data.match.ymax);
        mismatchScatterDetail = createChart("mismatchScatterDetail", data.mismatch.ymax);
        // mismatchScatterDetail = createChart("mismatchScatterDetail", data.mismatch.ymax);
    });
}

function afterSelection(event) {
    let label;
    event.preventDefault();
    var min = Math.trunc(event.xAxis[0].min);
    var max = Math.trunc(event.xAxis[0].max);

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

    coverageScatterMaster.xAxis[0].removePlotBand('mask-before');
    coverageScatterMaster.xAxis[0].addPlotBand({
        id: 'mask-before',
        from: min,
        to: max,
        color: 'rgba(0, 0, 0, 0.2)'
    });

    // // TODO this is where we will do the second AJAX call for selection
    // for (var i = 0; i < coverageScatterMaster.series[0].xData.length; i++) {
    //     let xAxisPoint = coverageScatterMaster.series[0].xData[i];
    //     if (xAxisPoint > min && xAxisPoint < max) {
    //         let yAxisPoint = coverageScatterMaster.series[0].yData[i];
    //         new_series.push([xAxisPoint, yAxisPoint]);
    //     }
    // }
    coverageScatterDetail.showLoading('Fetching data from the server');
    matchScatterDetail.showLoading('Fetching data from the server');
    $.getJSON('/api/v1/readuntil/benefitdata/detail', {min, max}, function (newSeriesDict, statusText, xhr) {
        console.log(newSeriesDict);
        if (xhr.status === 200) {
            coverageScatterDetail.series[0].setData(newSeriesDict.coverage.data);
            coverageScatterDetail.min = min;
            coverageScatterDetail.max = max;
            coverageScatterDetail.hideLoading();

            matchScatterDetail.series[0].setData(newSeriesDict.match.data);
            matchScatterDetail.min = min;
            matchScatterDetail.max = max;
            matchScatterDetail.hideLoading();
        } else {
            throw "Request error";
        }
    });


}

function drawReadUntilCharts() {
    $.getJSON('/api/v1/readuntil/benefitdata/master', function (data) {
        if (!Highcharts.Series.prototype.renderCanvas) {
            throw "Module not loaded";
        }

        let coverageDetail = $(".coverageScatterDetail").highcharts();
        // let coverageMaster = $(".coverageScatterMaster").highcharts();
        let min, max;
        // Wake up, Neo....
        // The matrix has you

        let axes;

        // Follow the white rabbit
        console.time("scatter");

        if (coverageDetail === undefined) {
            console.log("Creating for the first time");
            createCoverageMaster(data);
        } else {
            console.log("Redrawing for update");
            min = coverageDetail.min;
            max = coverageDetail.max;
            axes = [xaxisMin, xaxisMax];
            if min
            $.getJSON('/api/v1/readuntil/benefitdata/detail', {min, max}, function (newSeriesDict, statusText, xhr) {
                updateExistingChart(coverageScatterMaster, data);
                // TODO this is where we will do the second AJAX call for update
                updateExistingChart(coverageDetail, newSeriesDict.coverage.data);
                updateExistingChart(matchScatterDetail, newSeriesDict.match.data);
                updateExistingChart(mismatchScatterDetail, newSeriesDict.mismatch.data);
            });
        }
        console.timeEnd("scatter");


        // //////////////////////////////////////////////////
        // /////////// Local benefit scatter plot ///////////
        // //////////////////////////////////////////////////
        // Highcharts.chart("localBenScatterDetail", {
        //
        //     chart: {
        //         zoomType: "x",
        //         height: "400px"
        //     },
        //
        //     boost: {
        //         useGPUTranslations: true,
        //         usePreAllocated: true
        //     },
        //
        //     xAxis: {
        //         min: 0,
        //         max: data.localBenefit.xmax,
        //         gridLineWidth: 1
        //     },
        //
        //     yAxis: {
        //         // Renders faster when we don"t have to compute min and max
        //         min: 0,
        //         max: data.localBenefit.ymax,
        //         title: {
        //             text: null
        //         }
        //     },
        //
        //     title: {
        //         text: "Knock, knock, Neo."
        //     },
        //
        //     legend: {
        //         enabled: false
        //     },
        //
        //     series: [{
        //         type: "scatter",
        //         color: "orange",
        //         data: data.localBenefit.data,
        //         marker: {
        //             fillcolor: "red",
        //             radius: 1
        //         },
        //         tooltip: {
        //             followPointer: false
        //         }
        //     }]
        //
        // });
        // ////////////////////////////////////////
        // ///////// Rolling Benefit line /////////
        // ////////////////////////////////////////
        // Highcharts.chart("rollingBenefitLineDetail", {
        //
        //     chart: {
        //         zoomType: "x",
        //         height: "400px"
        //     },
        //
        //     boost: {
        //         useGPUTranslations: true,
        //         usePreAllocated: true
        //     },
        //
        //     xAxis: {
        //         min: 0,
        //         max: data.rollingBenefitFwd.xmax,
        //         gridLineWidth: 1
        //     },
        //
        //     yAxis: {
        //         // Renders faster when we don"t have to compute min and max
        //         min: 0,
        //         max: data.rollingBenefitFwd.ymax,
        //         title: {
        //             text: null
        //         }
        //     },
        //
        //     title: {
        //         text: "Neo: Who is it?"
        //     },
        //
        //     legend: {
        //         enabled: false
        //     },
        //
        //     series: [{
        //         type: "line",
        //         name: "Forward rolling benefit",
        //         color: "white",
        //         data: data.rollingBenefitFwd.data,
        //         marker: {
        //             fillcolor: "white",
        //             radius: 1
        //         }
        //     }, {
        //         type: "line",
        //         name: "Reverse rolling benefit",
        //         color: "gold",
        //         data: data.rollingBenefitRev.data,
        //         marker: {
        //             fillcolor: "gold",
        //             radius: 1
        //         }
        //         // },
        //         // tooltip: {
        //         //     followPointer: false,
        //         //     pointFormat: "[{point.x:.1f}, {point.y:.1f}]",
        //         //     valueDecimals: 8
        //         // }
        //     }]
        //
        // });
        // ///////////////////////////////////////////////
        // /////////// Mask Stepped Line Forward//////////
        // ///////////////////////////////////////////////
        // Highcharts.chart("maskStepLineFwdDetail", {
        //
        //     chart: {
        //         zoomType: "x",
        //         height: "400px"
        //     },
        //
        //     boost: {
        //         useGPUTranslations: true,
        //         usePreAllocated: true
        //     },
        //
        //     xAxis: {
        //         min: 0,
        //         max: data.forwardMask.xmax,
        //         gridLineWidth: 1
        //     },
        //
        //     yAxis: {
        //         // Renders faster when we don"t have to compute min and max
        //         min: 0,
        //         max: data.forwardMask.ymax,
        //         title: {
        //             text: null
        //         }
        //     },
        //
        //     title: {
        //         text: "Choi: It's Choi."
        //     },
        //
        //     legend: {
        //         enabled: false
        //     },
        //     series: [{
        //         type: "line",
        //         step: "right",
        //         name: "Forward mask",
        //         color: "gold",
        //         data: data.forwardMask.data,
        //         marker: {
        //             fillcolor: "gold",
        //             radius: 1
        //         }
        //     }]
        //
        // });
        // ////////////////////////////////////////////////
        // /////////// Mask Stepped Line Reverse //////////
        // ////////////////////////////////////////////////
        // Highcharts.chart("maskStepLineRevDetail", {
        //
        //     chart: {
        //         zoomType: "x",
        //         height: "400px"
        //     },
        //
        //     boost: {
        //         useGPUTranslations: true,
        //         usePreAllocated: true
        //     },
        //
        //     xAxis: {
        //         min: 0,
        //         max: data.forwardMask.xmax,
        //         gridLineWidth: 1
        //     },
        //
        //     yAxis: {
        //         // Renders faster when we don"t have to compute min and max
        //         min: 0,
        //         max: data.forwardMask.ymax,
        //         title: {
        //             text: null
        //         }
        //     },
        //
        //     title: {
        //         text: "Neo: You're two hours late."
        //     },
        //
        //     legend: {
        //         enabled: false
        //     },
        //     series: [{
        //         type: "line",
        //         name: "Reverse mask",
        //         color: "white",
        //         data: data.reverseMask.data,
        //         step: "right",
        //         marker: {
        //             fillcolor: "white",
        //             radius: 1
        //         }
        //     }]
        //
        // });

    });
}