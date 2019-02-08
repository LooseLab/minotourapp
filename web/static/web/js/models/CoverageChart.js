class CoverageChart {

    constructor(div_name) {

        this._div_master_name = div_name + "_master";
        this._div_detail_name = div_name + "_detail";

        this._master_chart = this.create_master_chart();
        this._detail_chart = this.create_detail_chart();
    }

    get master_chart() {

        return this._master_chart;
    }

    get detail_chart() {

        return this._detail_chart;
    }

    create_master_chart() {

        let self = this;

        var chart = new Highcharts.chart({

            chart: {
                renderTo: this._div_master_name,
                reflow: true,
                backgroundColor: null,
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

        return chart;
    }

    create_detail_chart() {

        let self = this;

        var chart = Highcharts.chart({

            chart: {
                renderTo: this._div_detail_name,
                reflow: false,
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
                enabled: false
            }
        });

        return chart;
    }

    afterSelection(event) {

        event.preventDefault();

        var min = Math.trunc(event.xAxis[0].min);
        var max = Math.trunc(event.xAxis[0].max);

        this._master_chart.xAxis[0].removePlotBand('mask-before');
        this._master_chart.xAxis[0].addPlotBand({
            id: 'mask-before',
            from: min,
            to: max,
            color: 'rgba(0, 0, 0, 0.2)'
        });

        var new_serie = [];

        for (var i = 0; i < this._master_chart.series[0].xData.length; i++) {

            if (this._master_chart.series[0].xData[i] > min && this._master_chart.series[0].xData[i] < max) {
                new_serie.push([this._master_chart.series[0].xData[i], this._master_chart.series[0].yData[i]]);
            }
        }

        this._detail_chart.series[0].setData(new_serie);
        this._detail_chart.min = min;
        this._detail_chart.max = max;
    }
}