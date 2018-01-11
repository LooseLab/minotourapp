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

function dynamicSort(property) {
    var sortOrder = 1;
    if (property[0] === "-") {
        sortOrder = -1;
        property = property.substr(1);
    }
    return function (a, b) {
        var result = (a[property] < b[property]) ? -1 : (a[property] > b[property]) ? 1 : 0;
        return result * sortOrder;
    }
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


function drawtaskbutton(taskstring, colour, icon, description, message, percentage, message2, i, long_description, reference, name, transcriptome) {
    taskstring = taskstring + '<div class="col-md-4">';
    taskstring = taskstring + '<button type="button" class="info-box ' + colour + '" data-toggle="modal" data-target="#taskmodal' + i + '">';
    taskstring = taskstring + '<div >';
    taskstring = taskstring + '<span class="info-box-icon"><i class="' + icon + '"></i></span>';
    taskstring = taskstring + '<div class="info-box-content">';
    taskstring = taskstring + '<span class="info-box-text">' + description + '</span>';

    taskstring = taskstring + '<span class="info-box-number" id="' + name + '-message">' + message + '</span>';
    taskstring = taskstring + '<div class="progress">';
    taskstring = taskstring + '<div class="progress-bar" id="' + name + '-percentage" style="width: ' + percentage + '%"></div>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '<span class="progress-description" id="' + name + '-message2" >';
    taskstring = taskstring + message2;
    taskstring = taskstring + '</span>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '</button>';
    taskstring = taskstring + '<div id="taskmodal' + i + '" class="modal fade" role="dialog">';
    taskstring = taskstring + '<div class="modal-dialog">';
    taskstring = taskstring + '<div class="modal-content">';
    taskstring = taskstring + '<div class="modal-header">';
    taskstring = taskstring + '<button type="button" class="close" data-dismiss="modal">&times;</button>';
    taskstring = taskstring + '<h4 class="modal-title">' + description + '</h4>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '<div class="modal-body">';
    taskstring = taskstring + '<p>' + long_description + '</p>';
    if (reference == true) {
        taskstring = taskstring + "<p>Select Reference:</p>"
        taskstring = taskstring + '<select id="' + name + '">';
        //console.log(this.references);
        for (var j = 0; j < this.references.length; j++) {
            if (this.references[j]['transcripts'] == transcriptome) {
                taskstring = taskstring + '<option>' + this.references[j]["reference_name"] + '</option>';
            }
        }
        taskstring = taskstring + '</select>';
    }
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '<div class="modal-footer">';
    taskstring = taskstring + '<button type="button" class="btn btn-default" id="button' + name + '">Go</button>';
    taskstring = taskstring + '<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '</div>';
    taskstring = taskstring + '</div>';

    return taskstring;
};

function write_run_data(textinfo) {
    console.log("textshiz");
    console.log(textinfo);
    var text = '';
    text += "<div class='table-responsive table-bordered'>";
    text += "<table class='table'>";
    text += "<tr>";
    text += "<th>Run</th>";
    text += "<th>Run Start Time</th>";
    text += "<th>Last Read Seen</th>";
    text += "<th>MinKNOW Computer Name</th>";
    text += "<th>MinION ID</th>";
    text += "<th>ASIC ID</th>";
    text += "<th>Current Script</th>";
    text += "<th>minKNOW version</th>";
    text += "<th>FlowCell Name</th> ";
    text += "<th>Sample Name</th>";
    text += "<th>Run Name</th>";
    text += "<th>FLow Cell ID</th>";
    text += "<th>Sequencing</th>";
    text += "<th>Barcoded</th>";


    text += "</tr>";
    $.each(textinfo, function (key, val) {
        text += "<tr>";
        text += "<td>";
        text += key;
        text += textinfo[key];
        text += "</td>";


        text += "<td>" + textinfo[key]['starttime'] + "</td>";
        text += "<td>" + textinfo[key]['lastread'] + "</td>";
        text += "<td>" + textinfo[key]["computer_name"] + "</td>";
        text += "<td>" + textinfo[key]["minIONname"] + "</td>";
        text += "<td>" + textinfo[key]["minKNOW_asic_id"] + "</td>";
        text += "<td>" + textinfo[key]["minKNOW_current_script"] + "</td>";
        text += "<td>" + textinfo[key]["minKNOW_version"] + "</td>";
        text += "<td>" + textinfo[key]['name'] + "</td>";
        text += "<td>" + textinfo[key]["sample_name"] + "</td>";

        text += "<td>" + textinfo[key]['run_name'] + "</td>";
        text += "<td>" + textinfo[key]['flowcellid'] + "</td>";
        if (textinfo[key]['active'] == true) {
            text += '<td><i class="fa fa-check" aria-hidden="true"></i></td>';
        } else {
            text += '<td><i class="fa fa-times" aria-hidden="true"></i></td>';
        }
        if (textinfo[key]['barcodes'].length > 2) {
            text += '<td><i class="fa fa-check" aria-hidden="true"></i></td>';
        } else {
            text += '<td><i class="fa fa-times" aria-hidden="true"></i></td>';
        }

        text += "</tr>";
    })
    text += "</table></div>";
    //console.log(text);
    $("#target_for_data").html(text);
    //$("#target_for_data").html(text);
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
            enabled: true
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
            min: 0,
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