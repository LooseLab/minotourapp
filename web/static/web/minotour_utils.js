Highcharts.setOptions({
    plotOptions: {
        series: {
            animation: false,
            turboThreshold: 0
        }
    }
});

/**
 * Get the selected flowcell tab
 */
function getSelectedTab() {
    return sessionStorage.getItem("flowcellTab");
}

/**
 * Set the selected flowcell Tab
 * @param tab {string} The tab that a user has just switched to
 */
function setSelectedTab(tab) {
    sessionStorage.setItem("flowcellTab", tab);
}

/**
 * Set the Flowcell ID in the session storage for whatever flowcell we are inspecting.
 * @param flowcell_id
 */
function setSelectedFlowcell(flowcellId) {
    sessionStorage.setItem("flowcellId", flowcellId);
}

/**
 * Return the selected Flowcell from the hidden input.
 */
function getSelectedFlowcell() {
    return sessionStorage.getItem("flowcellId");
}


/**
 * Return the barcode name for given tab from session storage.
 * @param tab {string} Tab name
 * @returns {string} Barcode name for this tab.
 */
function getSelectedBarcode(tab) {
    return sessionStorage.getItem(`${tab}-barcode`);
}

/**
 * Set the barcode of the tab the user is on in the session storage.
 * @param barcode {string} Barcode name
 * @param tab {string} Tab name
 */
function setSelectedBarcode(barcode, tab) {
    sessionStorage.setItem(`${tab}-barcode`, barcode);
}


/**
 * Remove old series from a given chart
 * @param chart {object}
 */
function clearChartData(chart) {
    while (chart.series.length > 0) {
        chart.series[0].remove();
    }
}


/**
 * Check whether the data we are about to add to a chart series is identical to the data already in that series.
 * @param newChartData {array} The data we plan to insert.
 * @param oldChartData {array} The data that exists in the series. Accessed via series.options.data.
 * @returns {bool} Whether the arrays are identical or not.
 */
function checkHighChartsDataIsNew(newChartData, oldChartData) {
    let flattenedNewData = [];
    let flattenedOldData = [];
    let identical = false;
    flattenedNewData = [].concat(...newChartData);
    // flatten chart data so we can compare it to the new data, can't compare nested arrays
    flattenedOldData = [].concat(...oldChartData);

    // compare old data with new data
    identical = flattenedNewData.length === flattenedOldData.length && flattenedNewData.every((value, index) => value === flattenedOldData[index]);
    return identical;
}

/**
 * Get a cookie value.
 * @param name {string} The name of the cookie to be parsed.
 * @returns {string} The value of the given cookie.
 */
function getCookie(name) {
    var cookieValue = null;
    if (document.cookie && document.cookie !== '') {
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = jQuery.trim(cookies[i]);
            // Does this cookie string begin with the name we want?
            if (cookie.substring(0, name.length + 1) === (name + '=')) {
                cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                break;
            }
        }
    }
    return cookieValue;
}


/**
 * @param {number} item_index The index of the navbar element in the
 * @function {function} Set the active item on the top navbar.
 *
 *
 */
function set_active_navbar_item(item_index) {

    var nav_bar = document.querySelectorAll('.navbar li a');

    nav_bar.forEach(function (element, index) {

        element.classList.remove("active");
        if (index == item_index) element.classList.add('active');
    });
}


/**
 * @function makeColumnChart
 * @param {string} divId ID of the div we want to draw the charts in
 * @param {string} chartTitle
 * @param {string} yAxisTitle
 * Creates an empty highCharts column chart at the provided divName,
 * with divName being the Id if the div we want to draw in.
 */
function makeColumnChart(divId, chartTitle, yAxisTitle) {

    var chart;
    chart = Highcharts.chart(divId, {
        chart: {
            type: "column",
            animation: Highcharts.svg, // don"t animate in old IE
            marginRight: 10,
            zoomType: "xy"
        },
        title: {
            text: chartTitle
        },
        yAxis: {
            title: {
                text: yAxisTitle
            },
            plotLines: [{
                value: 0,
                width: 1,
                color: "#808080"
            }]
        },
        lang: {
            noData: "Looking for data probably."
        },
        legend: {
            enabled: true
        },
        exporting: {
            enabled: true
        },
        series: []
    });

    return chart;
}


/**
 * Function abstracting creation of a HighCharts Spline chart.
 * @param divName {str} ID of div to initialise chart inside.
 * @param chartTitle {str} Title of chart.
 * @param yAxisTitle {str} Title of Y axis.
 * @returns {*}
 */
function makeSplineChart(divId, chartTitle, yAxisTitle) {
    const chart = Highcharts.chart(divId, {
        chart: {
            type: "spline",
            marginRight: 10,
            animation: false,
            zoomType: "x"
        },
        boost: {
            useGPUTranslations: true
        },
        title: {
            text: chartTitle
        },
        xAxis: {
            type: "datetime",
            tickPixelInterval: 150
        },
        yAxis: {
            title: {
                text: yAxisTitle
            },
            plotLines: [{
                value: 0,
                width: 1,
                color: "#808080"
            }]
        },
        legend: {
            enabled: true
        },
        exporting: {
            enabled: true,
            sourceWidth: 1200,
            sourceHeight: 400,
        }
    });

    return chart;
}


function makeSplineChartNonDatetime(divName, chartTitle, yAxisTitle, xAxisTitle) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: "spline",
            marginRight: 10,
            animation: false,
            zoomType: "x"
        },
        boost: {
            useGPUTranslations: true
        },
        title: {
            text: chartTitle
        },
        xAxis: {
            title: {
                text: xAxisTitle
            }
        },
        yAxis: {
            title: {
                text: yAxisTitle
            },
            plotLines: [{
                value: 0,
                width: 1,
                color: "#808080"
            }]
        },
        legend: {
            enabled: true
        },
        exporting: {
            enabled: true
        }
    });

    return chart;
}

function makeBoxPlot(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: 'boxplot'
        },
        title: {
            text: chartTitle
        },
        legend: {
            enabled: true
        },
        xAxis: {
            type: 'category'
        },
        yAxis: {
            title: {
                text: yAxisTitle
            },
            min: 0
        }
    });

    return chart;
}

function makeLiveHistogram(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: 'column',
            //marginRight: 10,
            animation: false,
            zoomType: 'x'
        },
        title: {
            text: chartTitle,
        },
        xAxis: {
            categories: [],
        },
        yAxis: {
            title: {
                text: 'Total Event Length'
            }
        },
        credits: {
            enabled: false
        },
        series: [{
            name: 'Read Histogram',
            //data: this.datain
        }]
    });
    return chart;

}

function makeAreaPlot(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.stockChart(divName, {
        chart: {
            renderTo: 'container-porehist' + this.title,
            type: 'area',
            //type: 'spline',
            height: 350,
            marginRight: 200,
        },
        title: {
            text: chartTitle,
        },
        xAxis: {
            range: 1 * 360 * 1000, //set range to last hour of data
        },
        rangeSelector: {
            enabled: false
        },
        yAxis: {
            //max: 512,
            endOnTick: false,
            title: {
                text: 'Channel Classifications'
            }
        },
        legend: {
            enabled: true,
            align: 'right',
            verticalAlign: 'top',
            layout: 'vertical',
            x: 0,
            y: 50
        },
        plotOptions: {
            area: {
                stacking: 'percent',
            },
            series: {
                showInNavigator: true,
                dataLabels: {
                    enabled: false,
                    formatter: function () {
                        return this.y;
                    }
                }
            }
        },
        credits: {
            enabled: false
        },
        series: [],
    });
    return chart;
}


function makeLiveChart(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.stockChart(divName, {
        chart: {
            type: 'spline',
            zoomType: 'xy'
        },
        boost: {
            useGPUTranslations: true
        },
        rangeSelector: {
            enabled: true,
            buttons: [{
                type: 'minute',
                count: 1,
                text: '1min'
            }, {
                type: 'minute',
                count: 5,
                text: '5min'
            }, {
                type: 'minute',
                count: 30,
                text: '1/2hr'
            }, {
                type: 'minute',
                count: 60,
                text: '1hr'
            }, {
                type: 'day',
                count: 0.5,
                text: '12hrs'
            }, {
                type: 'day',
                count: 1,
                text: '1day'
            }, {
                type: 'all',
                text: 'All'
            }]
        }
        ,
        title: {
            text: chartTitle,
        }
        ,
        xAxis: {
            type: 'datetime',
            tickPixelInterval: 150
        }
        ,
        yAxis: {
            title: {
                text: yAxisTitle
            }
            ,
            plotLines: [{
                value: 0,
                width: 1,
                color: '#808080'
            }],
            //min: 0,
        }
        ,
        credits: {
            enabled: false
        }
        ,
        series: []
    });
    return chart;
};

/**
 * Create a HighCharts heatmap chart. These are used on the basecalled data tab.
 * @param divName {string} The id of the div to render chart to.
 * @param chartTitle {string} The title of the chart.
 * @return {*}
 */
function makeHeatmapChart(divName, chartTitle) {

    var chart = Highcharts.chart(divName, {
        chart: {
            type: "heatmap",
        },
        title: {
            text: chartTitle
        },
        colorAxis: {
            min: 0,
            minColor: '#FFFFFF',
            maxColor: Highcharts.getOptions().colors[0]
        },
        xAxis: {
            title: null,
            labels: {
                enabled: false
            },

        },
        yAxis: {
            title: null,
            labels: {
                enabled: false
            },

        },
        legend: {
            layout: 'horizontal',
        },
        exporting: {
            enabled: false
        }
    });

    return chart;

}

/**
 * Create a HighCharts pie chart. These are used on the basecalled data tab.
 * @param divId {string} The id of the div to render chart to.
 * @param chartTitle {string} The title of the chart.
 * @return {*}
 */
function makePieChart(divId, chartTitle) {
    let chart = Highcharts.chart(divId, {
        chart: {
            plotBackgroundColor: null,
            plotBorderWidth: null,
            plotShadow: false,
            type: 'pie'
        },
        title: {
            text: chartTitle
        },
        tooltip: {
            pointFormat: '{series.name}: <b>{point.percentage:.1f}%</b>'
        },
        accessibility: {
            point: {
                valueSuffix: '%'
            }
        },
        lang: {
            noData: "Looking for data probably."
        },
        plotOptions: {
            pie: {
                allowPointSelect: true,
                cursor: 'pointer',
                dataLabels: {
                    enabled: true,
                    format: '<b>{point.name}</b>: {point.percentage:.1f} %'
                }
            }
        },
        series: []
    })
    return chart
}
