Highcharts.setOptions({
    plotOptions: {
        series: {
            animation: false,
            turboThreshold: 0
        }
    }
});

var NUMBER_SECONDS_IN_A_MINUTE = 60;

var MINOTOUR_VERSION = 0.5;

function csrfSafeMethod(method) {
    return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
}

function check_minotour_version() {
    $.getJSON("http://www.nottingham.ac.uk/~plzloose/minoTourhome/message.php?callback=?", function (result) {

        $.each(result, function (key, value) {
            //checking version info.
            if (key == "version") {
                if (value == MINOTOUR_VERSION) {
                    $("#newstarget").html("You are running the most recent version of minoTour - version " + value + ".<br>");
                }
                else if (value < MINOTOUR_VERSION) {
                    $("#newstarget").html("You appear to be in the fortunate position of running a future version of the minoTour web application " + value + ". If you have modified the code yourself - great. If not then there might be an issue somewhere!.<br>");
                }
                else if (value > MINOTOUR_VERSION) {
                    $("#newstarget").html("You are running an outdated version of the minoTour web application. The most recent version of minoTour is version " + value + ".<br>" + "Instructions for upgrading will be posted below.<br>");
                }

            } else if (key.substring(0, 7) == "message") {
                $("#newstarget").append(value + "<br>");
            }
        });
    });
}


function check_user_runs() {
    var url = "/api/v1/runs/";

    $.getJSON(url, function (data) {
        var items = [];
        $.each(data, function (key, val) {
            items.push("<li id=" + key + "><a href='/web/private/runs/" + val.id + "'>" + val.run_name + "</a></li>");
        });

        var text = "<span>You have " + data.length + " minION runs available to view.</span>";
        text += "<ul>";
        text += items.join("");
        text += "</ul>";

        $("#div-user-runs").html(text);
    });
}


function makeChart(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: "column",
            animation: Highcharts.svg, // don"t animate in old IE
            marginRight: 10
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

function makeChartlabels(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: "column",
            animation: Highcharts.svg, // don"t animate in old IE
            marginRight: 10,
            zoomType: "x"
        },
        title: {
            text: chartTitle
        },
        xAxis: {
            type: 'category',
            labels: {
                rotation: -90,
                style: {
                    fontSize: '9px',
                    fontFamily: 'Verdana, sans-serif'
                }
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
        },
        series: []
    });

    return chart;
}

function makeChart2(divName, chartTitle, yAxisTitle) {
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


function makeChart3(divName, chartTitle, yAxisTitle, xAxisTitle) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: "scatter",
            marginRight: 10,
            animation: false,
            zoomType: "xy"
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

function makeChart4(divName, chartTitle, yAxisTitle, xAxisTitle) {
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

function makeYieldProjection(divName, chartTitle, yAxisTitle) {
    var chart = Highcharts.chart(divName, {
        chart: {
            type: 'spline',
            zoomType: 'x',
            height: 350,
        },
        boost: {
            useGPUTranslations: true
        },
        title: {
            text: chartTitle,
        },
        xAxis: {
            type: 'datetime',
            tickPixelInterval: 150
        },
        yAxis: {
            title: {
                text: 'Cumulative Predicted Bases'
            },
            plotLines: [{
                value: 0,
                width: 1,
                color: '#808080'
            }],
            min: 0,
        },
        credits: {
            enabled: false
        },
        series: [{
            name: 'Original Data',
            data: []
        },
            {
                name: 'Current Projected Data',
                dashStyle: 'longdash',
                data: []
            },
            {
                name: 'First Hour Projected Data',
                dashStyle: 'shortdash',
                data: []
            },
            {
                name: 'Ideal Results',
                dashStyle: 'Dot',
                data: []
            }
        ]

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

//
// https://stackoverflow.com/questions/10343838/get-value-of-span-text
//
function copyToClipboardMsg(elem, msgElem) {
	  var succeed = copyToClipboard(elem);
    var msg;
    if (!succeed) {
        msg = "Copy not supported or blocked.  Press Ctrl+c to copy."
    } else {
        msg = "Text copied to the clipboard."
    }
    if (typeof msgElem === "string") {
        msgElem = document.getElementById(msgElem);
    }
    msgElem.innerHTML = msg;
    setTimeout(function() {
        msgElem.innerHTML = "";
    }, 3000);
}

function copyToClipboard(elem) {
	  // create hidden text element, if it doesn't already exist
    var targetId = "_hiddenCopyText_";
    var isInput = elem.tagName === "INPUT" || elem.tagName === "TEXTAREA";
    var origSelectionStart, origSelectionEnd;
    if (isInput) {
        // can just use the original source element for the selection and copy
        target = elem;
        origSelectionStart = elem.selectionStart;
        origSelectionEnd = elem.selectionEnd;
    } else {
        // must use a temporary form element for the selection and copy
        target = document.getElementById(targetId);
        if (!target) {
            var target = document.createElement("textarea");
            target.style.position = "absolute";
            target.style.left = "-9999px";
            target.style.top = "0";
            target.id = targetId;
            document.body.appendChild(target);
        }
        target.textContent = elem.textContent;
    }
    // select the content
    var currentFocus = document.activeElement;
    target.focus();
    target.setSelectionRange(0, target.value.length);

    // copy the selection
    var succeed;
    try {
    	  succeed = document.execCommand("copy");
    } catch(e) {
        succeed = false;
    }
    // restore original focus
    if (currentFocus && typeof currentFocus.focus === "function") {
        currentFocus.focus();
    }

    if (isInput) {
        // restore prior selection
        elem.setSelectionRange(origSelectionStart, origSelectionEnd);
    } else {
        // clear temporary content
        target.textContent = "";
    }
    return succeed;
}
