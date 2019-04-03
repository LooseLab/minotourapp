function drawCoverageScatter() {
    $.getJSON('/api/v1/readuntil/benefitdata?', function (data) {

        if (!Highcharts.Series.prototype.renderCanvas) {
            throw "Module not loaded";
        }
        // Wake up, Neo....
        // The matrix has you
        let coverageScatterDetail;
        let coverageScatterMaster;
        // Follow the white rabbit
        console.time("scatter");

        function createCoverageDetail(masterChart) {

            var detailData = [],
                detailStart = data.coverage.data[0][0];
            // console.log(masterChart.series[0]);
            // for (let i = 0; i < masterChart.series[0].xData.length; i++) {
            //     let xAxisPoint = masterChart.series[0].xData[i];
            //     let yAxisPoint = masterChart.series[0].yData[i];
            //     if (xAxisPoint >= detailStart) {
            //         detailData.push([xAxisPoint, yAxisPoint]);
            //     }
            // }
            coverageScatterDetail = Highcharts.chart("coverageScatterDetail", {

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
                    min: 0,
                    max: data.coverage.xmax,
                    gridLineWidth: 1
                },

                yAxis: {
                    // Renders faster when we don"t have to compute min and max
                    min: 0,
                    max: data.coverage.ymax,
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
        }

        function createCoverageMaster() {
            coverageScatterMaster = Highcharts.chart("coverageScatterMaster", {

                chart: {
                    zoomType: "x",
                    height: "150px",
                    marginLeft: 50,
                    marginRight: 20,
                    events: {
                        selection: afterSelection.bind(this)
                        // listen to the selection event on the master chart to update the
                        // extremes of the detail chart
                        // selection: function (event) {
                        //     var extremesObject = event.xAxis[0],
                        //         min = Math.trunc(extremesObject.min),
                        //         max = Math.trunc(extremesObject.max),
                        //         detailData = [],
                        //         xAxis = this.xAxis[0];
                        //     console.log(extremesObject);
                        //
                        //     // reverse engineer the last part of the data
                        //     console.log(this.series);
                        //     for (let i = 0; i < this.series[0].xData.length; i++) {
                        //         let xAxisPoint = this.series[0].xData[i];
                        //         let yAxisPoint = this.series[0].yData[i];
                        //         if (xAxisPoint > min && xAxisPoint < max) {
                        //             detailData.push([xAxisPoint, yAxisPoint]);
                        //         }
                        //     }
                        //     // move the plot bands to reflect the new detail span
                        //     xAxis.removePlotBand('mask-before');
                        //     xAxis.addPlotBand({
                        //         id: 'mask-before',
                        //         from: data.coverage.data[0][0],
                        //         to: min,
                        //         color: 'rgba(0, 0, 0, 0.2)'
                        //     });
                        //
                        //     xAxis.removePlotBand('mask-after');
                        //     xAxis.addPlotBand({
                        //         id: 'mask-after',
                        //         from: max,
                        //         to: data.coverage.data[data.coverage.data.length - 1][0],
                        //         color: 'rgba(0, 0, 0, 0.2)'
                        //     });
                        //
                        //
                        //     coverageScatterDetail.series[0].setData(detailData);
                        //     console.log(coverageScatterDetail);
                        //     console.log([max, min]);
                        //     coverageScatterDetail.min = min;
                        //     coverageScatterDetail.max = max;
                        //
                        //     return false;
                        // }
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

                title: {
                    text: "Wake up, Neo..."
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
                createCoverageDetail(masterChart);
            });
        }

        function updateExistingCoverageDetail(detailChart){

            detailChart.setData()
        }

        function updateExisitngCoverageMaster(masterChart){

        }
        function afterSelection(event) {

            event.preventDefault();

            var min = Math.trunc(event.xAxis[0].min);
            var max = Math.trunc(event.xAxis[0].max);

            coverageScatterMaster.xAxis[0].removePlotBand('mask-before');
            coverageScatterMaster.xAxis[0].addPlotBand({
                id: 'mask-before',
                from: min,
                to: max,
                color: 'rgba(0, 0, 0, 0.2)'
            });

            var new_series = [];

            for (var i = 0; i < coverageScatterMaster.series[0].xData.length; i++) {
                let xAxisPoint = coverageScatterMaster.series[0].xData[i];
                if (xAxisPoint > min && xAxisPoint < max) {
                    let yAxisPoint = coverageScatterMaster.series[0].yData[i];
                    new_series.push([xAxisPoint, yAxisPoint]);
                }
            }

            coverageScatterDetail.series[0].setData(new_series);
            coverageScatterDetail.min = min;
            coverageScatterDetail.max = max;
            console.log([min, max]);
            console.log(coverageScatterDetail.min);
            console.log(coverageScatterDetail.max);
        }

        createCoverageMaster();
        console.timeEnd("scatter");
        ///////////////////////////////////////////
        ///// Draw the matches scatter plot ///////
        ///////////////////////////////////////////
        // Highcharts.chart("matchScatterDetail", {
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
        //         max: data.match.xmax,
        //         gridLineWidth: 1
        //     },
        //
        //     yAxis: {
        //         // Renders faster when we don"t have to compute min and max
        //         min: 0,
        //         max: data.match.ymax,
        //         title: {
        //             text: null
        //         }
        //     },
        //
        //     title: {
        //         text: "The Matrix has you..."
        //     },
        //
        //     legend: {
        //         enabled: false
        //     },
        //
        //     series: [{
        //         type: "scatter",
        //         color: "green",
        //         data: data.match.data,
        //         marker: {
        //             fillcolor: "red",
        //             radius: 1
        //         },
        //         tooltip: {
        //             followPointer: false,
        //             pointFormat: "[{point.x:.1f}, {point.y:.1f}]",
        //             valueDecimals: 8
        //         }
        //     }]
        //
        // });
        // /////////////////////////////////////////////////
        // ////////// Draw the mismatch scatter ////////////
        // /////////////////////////////////////////////////
        // Highcharts.chart("mismatchScatterDetail", {
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
        //         max: data.mismatch.xmax,
        //         gridLineWidth: 1
        //     },
        //
        //     yAxis: {
        //         // Renders faster when we don"t have to compute min and max
        //         min: 0,
        //         max: data.mismatch.ymax,
        //         title: {
        //             text: null
        //         }
        //     },
        //
        //     title: {
        //         text: "Follow the white rabbit."
        //     },
        //
        //     legend: {
        //         enabled: false
        //     },
        //
        //     series: [{
        //         type: "scatter",
        //         color: "red",
        //         data: data.mismatch.data,
        //         marker: {
        //             fillcolor: "red",
        //             radius: 1
        //         },
        //         tooltip: {
        //             followPointer: false,
        //             pointFormat: "[{point.x:.1f}, {point.y:.1f}]",
        //             valueDecimals: 8
        //         }
        //     }]
        //
        // });
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