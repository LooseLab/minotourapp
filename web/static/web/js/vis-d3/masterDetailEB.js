// Initialise the plots so they are global variables
let coverageScatterDetail;
let coverageScatterMaster;
let costForwardScatterDetail;
let costReverseScatterDetail;
let fixedBenefitForwardDetail;
let fixedBenefitReverseDetail;
let scoresReverseDetail;
let scoresForwardDetail;
let forwardMaskDetail;
let revMaskDetail;
let benefitDetail;


// showWidth is the maximum number of bases we can select to view on the master chart
let showWidth = 250000;

let chromosome_chosen;

function draw_selected_chromsome(flowcellId, chromosome_chosen){
    /*
    Reload a new chromosomes data
    */
    $.getJSON('/api/v1/readuntil/benefitdata/master', {flowcellId, chromosome_chosen}, function (data) {
        coverageScatterMaster.yAxis[0].setExtremes(0, data.coverage.ymax);
        updateExistingChart(coverageScatterMaster, data.coverage.data, 0);

    });
}

function update_axes_and_data(chart, min, max, data, ymax){
    // Set the xAxis extremes to the newly selected values
    chart.xAxis[0].setExtremes(min, max);
    // if we need to update the yaxis as well
    if(ymax && data.ymax < 0){
        // if we have a ymax below 0 we need to flip the yaxis so negative at x axis and 0 at top
        console.log(chart);
        console.log(data.ymax-1);
        chart.yAxis[0].setExtremes(data.ymax, 0);
    }
    else if (ymax){
        // set the yAxisextremes
        chart.yAxis[0].setExtremes(0, data.ymax);
    }
    // set the data for the chart to display
    chart.series[0].setData(data.data);
    // update the axes
    chart.xAxis[0].update();

    // hide the loading screen
    chart.hideLoading();
}


function create_eb_selection_box(flowcellId) {
    // For all the sel objects that finding the element by the sel class returned
    $('.sel').each(function () {
        // get the child of the div, which is the select option in this case, do not display it
        $(this).children('select').css('display', 'none');
        // Get the sel div element as a var
        var $current = $(this);
        // get all the options children
        $.get("/api/v1/readuntil/benefitdata/chromosomes", {flowcellId}, function (chromosomes, statusText, xhr) {
            // if it's the first options clause
            if (xhr.status != 200) console.log("Error");
            // If we're on our first drawing of the select spans
            let first_iteration = !document.getElementById("eb-placeholder");
            // The chromosomes we have that are already present
            let alreadyPresentChromo = [];
            // If it is not our first drawing of it
            if (!first_iteration) {
                // The already present chromosomes html elements
                alreadyPresentChromo = [...$(".sel__box")[0].children];
                // for each of the elements
                alreadyPresentChromo.forEach(function(presChrom){
                    // if the chromsome in already present chromsome is in the list of chromosomes fetched from the server
                    if (chromosomes.includes(presChrom.textContent)){
                        // get the index,
                        let index = chromosomes.indexOf(presChrom.textContent);
                        // remove it from the freshly fetched chromosomes
                        chromosomes.splice(index, 1);
                    }
                });
                // If there are no new chromosomes
                if (!chromosomes.length){
                    return;
                }
            // if it is our first iteration
            } else {
                // unshift choose chromosomes
                chromosomes.unshift("Choose Chromosome");
            }


            chromosomes.forEach(function (chromosome, index) {

                // if this is the first time so we're creating the arrays
                if (index == 0 && first_iteration) {
                    // the div to the front of the sel div, make it sel__box class instead of sel
                    $current.prepend($('<div>', {
                        class: $current.attr('class').replace(/sel/g, 'sel__box ')

                    }));
                    // Set the placeholder text to the option text - This would be the first chromosome
                    var placeholder = chromosome;
                    // prepend a span in front of the box div, with the placeholder text
                    $current.prepend($('<span>', {
                        class: $current.attr('class').replace(/sel/g, 'sel__placeholder '),
                        text: placeholder,
                        id: "eb-placeholder",
                        'data-placeholder': placeholder
                    }));

                    return;
                }

                // Then for the rest of the options add them in
                $current.children('div').append($('<span>', {
                    class: $current.attr('class').replace(/sel/g, 'sel__box__options '),
                    text: chromosome
                }));
            });
            // Toggling the `.active` state on the `.sel`.
            $('.sel').click(function () {
                $(this).toggleClass('active');
            });

            // Toggling the `.selected` state on the options.
            $('.sel__box__options').click(function () {
                var txt = $(this).text();
                chromosome_chosen = txt;
                var index = $(this).index();


                $(this).siblings('.sel__box__options').removeClass('selected');
                $(this).addClass('selected');

                var $currentSel = $(this).closest('.sel');
                $currentSel.children('.sel__placeholder').text(txt);
                $currentSel.children('select').prop('selectedIndex', index + 1);

                draw_selected_chromsome(flowcellId, txt);


            });
            $(".sel").css("display", "");
        });
    });
}


function updateExistingChart(chart, newSeries, index) {
    // Set the data on an existing series on a chart, either master or detail
    chart.series[index].setData(newSeries);

}

////////////////////////////////////////////////////////////////////
///////// Draw the detail charts under the coverage detail /////////
////////////////////////////////////////////////////////////////////

function createChart(targetDivId, ymax, title, chartType = "scatter", stepped = false, legend = false) {
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

//////////////////////////////////////////////////////////////////////////
///////// Draw the single master chart under the coverage detail /////////
//////////////////////////////////////////////////////////////////////////
function createCoverageMaster(data, flowcellId) {
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
        coverageScatterDetail = createChart("coverageScatterDetail", data.coverage.ymax, "Coverage");
        benefitDetail = createChart("benefitDetail", 1, "Benefits");
        costForwardScatterDetail = createChart("costForwardScatterDetail", 1, "Forward cost", "line");
        costReverseScatterDetail = createChart("costReverseScatterDetail", 1, "Reverse cost", "line");
        scoresForwardDetail = createChart("scoresForwardDetail", 1, "Forward scores");
        scoresReverseDetail = createChart("scoresReverseDetail", 1, "Reverse scores");
        forwardMaskDetail = createChart("maskStepLineFwdDetail", 1, "Forward mask", "line", true);
        revMaskDetail = createChart("maskStepLineRevDetail", 1, "Reverse mask", "line", true);
    });
}

///////////////////////////////////////////////////////////////////////////////
///////// After we have selected an area on the master or detail chart ////////
///////////////////////////////////////////////////////////////////////////////

function afterSelection(event) {
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
    let charts = [coverageScatterDetail, benefitDetail, costForwardScatterDetail,
        costReverseScatterDetail,forwardMaskDetail, forwardMaskDetail,
        scoresForwardDetail, scoresReverseDetail];

    charts.forEach(function (chart) {
        chart.showLoading('Fetching data from the server');
    });

    // Get the detail chart data for the zoom selection from the server
    $.getJSON('/api/v1/readuntil/benefitdata/detail', {
        min,
        max,
        chromosome_chosen,
        flowcellId
    }, function (newSeriesDict, statusText, xhr) {

        // if successful request, set the detail charts to show the newly retrieve data
        if (xhr.status === 200) {

            update_axes_and_data(coverageScatterDetail, min, max, newSeriesDict.coverage, true);

            update_axes_and_data(benefitDetail, min, max, newSeriesDict.benefits, true);

            update_axes_and_data(forwardMaskDetail, min, max, newSeriesDict.forwardMask, false);

            update_axes_and_data(revMaskDetail, min, max, newSeriesDict.reverseMask, false);

            update_axes_and_data(costForwardScatterDetail, min, max, newSeriesDict.costsFwd, true);

            update_axes_and_data(costReverseScatterDetail, min, max, newSeriesDict.costsRev, true);

            update_axes_and_data(scoresForwardDetail, min, max, newSeriesDict.scoresFwd, true);

            update_axes_and_data(scoresReverseDetail, min, max, newSeriesDict.scoresRev, true);

        } else {
            throw "Request error";
        }
    });


}

// draw the charts, this is the function that is called in request data every 60 seconds
function drawReadUntilCharts() {
    let flowcellId = document.querySelector("#flowcell-id").value;
    create_eb_selection_box(flowcellId);

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

    let startData = {"coverage": {"ymax": 1, "data":[0,1]}};
    // If the coverage detail chart doesn't exist, draw it for the first time by calling the create master chart which has a callback that draws the other charts
    if (coverageDetail === undefined) {
        createCoverageMaster(startData, flowcellId);
    } else {
        // If the coverage detail chart does exist, update the existing charts
        // Get the minimum and maximum values of the xAxis
        min = coverageDetail.min;
        max = coverageDetail.max;
        // If min is not undefined, the charts are already displaying data, so update thew, otherwise we're good
        if (min !== undefined) {
            $.getJSON('/api/v1/readuntil/benefitdata/detail', {
                min,
                max,
                chromosome_chosen,
                flowcellId
            }, function (newSeriesDict, statusText, xhr) {
                updateExistingChart(coverageScatterMaster, data.coverage.data, 0);
                updateExistingChart(coverageDetail, newSeriesDict.coverage.data, 0);
                updateExistingChart(costForwardScatterDetail, newSeriesDict.costsFwd.data, 0);
                updateExistingChart(costReverseScatterDetail, newSeriesDict.costsRev.data, 0);
                updateExistingChart(scoresForwardDetail, newSeriesDict.scoresFwd.data, 0);
                updateExistingChart(scoresReverseDetail, newSeriesDict.scoresRev.data, 0);
                // updateExistingChart(costScatterDetail, newSeriesDict.match.data, 0);
                // updateExistingChart(mismatchScatterDetail, newSeriesDict.mismatch.data, 0);
                // updateExistingChart(localBenefitDetail, newSeriesDict.localBenefit.data, 0);
                // updateExistingChart(rollingBenefitLineDetail, newSeriesDict.rollingBenefitFwd.data, 0);
                // updateExistingChart(rollingBenefitLineDetail, newSeriesDict.rollingBenefitRev.data, 1);
                updateExistingChart(forwardMaskDetail, newSeriesDict.forwardMask.data, 0);
                updateExistingChart(revMaskDetail, newSeriesDict.reverseMask.data, 0);
            });
        }


    }
    console.timeEnd("scatter");

}