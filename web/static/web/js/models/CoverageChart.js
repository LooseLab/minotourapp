class CoverageChart {
    /*
    Create a new coverage chart, which is a highcharts stepped line chart
    returns a getter, and defines the master chart detail chart after selection event
     */

    constructor(div_name) {
        // the chart
        this._divIdMaster = div_name + "_master";
        this._divIdDetail = div_name + "_detail";

        this._masterChart = this.createMasterChart();
        this._detailChart = this.createDetailChart();
    }

    get masterChart() {

        return this._masterChart;
    }

    get detailChart() {

        return this._detailChart;
    }

    createMasterChart() {

        let self = this;

        return Highcharts.chart({

            chart: {
                renderTo: this._divIdMaster,
                reflow: true,
                marginLeft: 50,
                marginRight: 20,
                zoomType: 'x',
                step: true,
                events: {
                    selection: self.afterSelection.bind(self)
                }
            },
            title: {
                text: null
            },
            xAxis: {
                type: 'line',
                showLastTickLabel: true,
                plotBands: [{
                    id: 'mask-before',
                    color: 'rgba(0, 0, 0, 0.2)'
                }],
                title: {
                    text: null
                }
            },
            yAxis: {
                title: {
                    text: null
                },
                visible: false
            },
            legend: {
                enabled: false
            },
            credits: {
                enabled: false
            },
            plotOptions: {
                series: {
                    fillColor: {
                        // linearGradient: [0, 0, 0, 70],
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
                type: 'area',
                name: null,
                data: [],
                step: true
            }],
            exporting: {
                enabled: false
            }
        });
    }

    createDetailChart() {
        let self = this;
        return Highcharts.chart({
            chart: {
                renderTo: this._divIdDetail,
                reflow: true,
                marginLeft: 50,
                marginRight: 20,
                style: {
                    position: 'absolute'
                },
                step: true,
                zoomType: 'x',
                events: {
                    selection: self.afterSelection.bind(self)
                },
                type: 'area'
            },
            credits: {
                enabled: false
            },
            title: {
                text: null
            },
            xAxis: {
                type: 'line',
            },
            yAxis: {
                title: {
                    text: null
                },
                maxZoom: 0.1,
                min: 0
            },
            legend: {
                enabled: false
            },
            plotOptions: {
                series: {
                    marker: {
                        enabled: false,
                        states: {
                            hover: {
                                enabled: true,
                                radius: 3
                            }
                        }
                    }
                }
            },
            series: [{
                name: null,
                data: [],
                step: true
            }],
            exporting: {
                enabled: true
            }
        });

    }

    afterSelection(event) {
        let min, max, newDetailSeries;
        event.preventDefault();
        min = Math.trunc(event.xAxis[0].min);
        max = Math.trunc(event.xAxis[0].max);
        this._masterChart.xAxis[0].removePlotBand('mask-before');
        this._masterChart.xAxis[0].addPlotBand({
            id: 'mask-before',
            from: min,
            to: max,
            color: 'rgba(0, 0, 0, 0.2)'
        });
        newDetailSeries = this._masterChart.series[0].options.data.slice(min, max + 1)
        for (let i = min; i < this._masterChart.series[0].xData.length; i++) {
            if (this._masterChart.series[0].xData[i] < max) {
                newDetailSeries.push([this._masterChart.series[0].options.data[i]]);
            }
        }
        this._detailChart.series[0].setData(newDetailSeries);
        this._detailChart.min = min;
        this._detailChart.max = max;
    }
}