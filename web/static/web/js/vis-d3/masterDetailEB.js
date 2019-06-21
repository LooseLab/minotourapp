// Initialise the plots so they are global variables
let coverageScatterDetail;
let coverageScatterMaster;
let matchScatterDetail;
let mismatchScatterDetail;
let localBenefitDetail;
let rollingBenefitLineDetail;
let forwardMaskDetail;
let revMaskDetail;
let benefitDetail;


// showWidth is the maximum number of bases we can select to view on the master chart
let showWidth = 250000;

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
            // add choose chromosome to start of array to be the placeholder

            let first_iteration = !document.getElementById("eb-placeholder");
            console.log("first iteration is ");
            console.log(first_iteration);
            let alreadyPresentChromo = [];

            if (!first_iteration) {
                alreadyPresentChromo = [...$(".sel__box")[0].children];
                alreadyPresentChromo.forEach(function(presChrom){
                    if (chromosomes.includes(presChrom.textContent)){
                        let index = chromosomes.indexOf(presChrom.textContent);
                        chromosomes.splice(index, 1);
                    };
                });

                if (!chromosomes.length){
                    console.log("No new chromosomes");
                    return;
                }
            } else {
                chromosomes.unshift("Choose Chromosome");
            }
                            console.log("chromosomes are ");
                console.log(chromosomes);


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
                var index = $(this).index();

                $(this).siblings('.sel__box__options').removeClass('selected');
                $(this).addClass('selected');

                var $currentSel = $(this).closest('.sel');
                $currentSel.children('.sel__placeholder').text(txt);
                $currentSel.children('select').prop('selectedIndex', index + 1);


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
        // matchScatterDetail = createChart("matchScatterDetail", data.match.ymax, "Matches");
        // mismatchScatterDetail = createChart("mismatchScatterDetail", data.mismatch.ymax, "Mismatches");
        // localBenefitDetail = createChart("localBenScatterDetail", data.localBen.ymax, "Local benefit");
        // rollingBenefitLineDetail = createChart("rollingBenefitLineDetail", data.forwardRoll.ymax, "Rolling benefit", "line", false, true);
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
    coverageScatterDetail.showLoading('Fetching data from the server');
    benefitDetail.showLoading('Fetching data from the server');
    // matchScatterDetail.showLoading('Fetching data from the server');
    // mismatchScatterDetail.showLoading('Fetching data from the server');
    // localBenefitDetail.showLoading('Fetching data from the server');
    // rollingBenefitLineDetail.showLoading('Fetching data from the server');
    forwardMaskDetail.showLoading('Fetching data from the server');
    revMaskDetail.showLoading('Fetching data from the server');

    // Get the detail chart data for the zoom selection from the server
    $.getJSON('/api/v1/readuntil/benefitdata/detail', {
        min,
        max,
        flowcellId
    }, function (newSeriesDict, statusText, xhr) {

        // if successful request, set the detail charts to show the newly retrieve data
        if (xhr.status === 200) {
            console.log(min, max);
            console.log(Number.isInteger(min));
            coverageScatterDetail.xAxis[0].setExtremes(min, max);
            // coverageScatterDetail.xAxis[0].min = min;
            coverageScatterDetail.series[0].setData(newSeriesDict.coverage.data);
            coverageScatterDetail.xAxis[0].max = max;
            coverageScatterDetail.hideLoading();

            // check if we have a reverse line already or not on the rolling benefit detail chart
            // if (rollingBenefitLineDetail.series[1] === undefined){
            //     // if there is no line, add a new series to the forward reverse rolling benefit line for the reverse line data
            //     rollingBenefitLineDetail.addSeries({name: "Reverse rolling benefit",
            //     data:newSeriesDict.rollingBenefitRev.data, "color": "gold"});
            // } else {
            //     // if there is a line, update the already existing series on the forward line chart
            //     rollingBenefitLineDetail.series[1].setData(newSeriesDict.rollingBenefitRev.data);
            // }

            benefitDetail.xAxis[0].setExtremes(min, max);
            benefitDetail.series[0].setData(newSeriesDict.benefits.data);
            benefitDetail.hideLoading();


            forwardMaskDetail.xAxis[0].setExtremes(min, max);
            forwardMaskDetail.series[0].setData(newSeriesDict.forwardMask.data);
            forwardMaskDetail.hideLoading();

            revMaskDetail.series[0].setData(newSeriesDict.reverseMask.data);
            revMaskDetail.xAxis[0].setExtremes(min, max);
            revMaskDetail.hideLoading();

        } else {
            throw "Request error";
        }
    });


}

// draw the charts, this is the function that is called in request data every 60 seconds
function drawReadUntilCharts(chromsome) {
    let flowcellId = document.querySelector("#flowcell-id").value;
    create_eb_selection_box(flowcellId);
    $.getJSON('/api/v1/readuntil/benefitdata/master', {flowcellId}, function (data) {
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
            createCoverageMaster(data, flowcellId);
        } else {
            // If the coverage detail chart does exist, update the existing charts
            // Get the minimum and maximum values of the xAxis
            min = coverageDetail.min;
            max = coverageDetail.max;
            // If min is not undefined, the charts are already displaying data, so update the, otherwise we're good
            if (min !== undefined) {
                $.getJSON('/api/v1/readuntil/benefitdata/detail', {
                    min,
                    max,
                    flowcellId
                }, function (newSeriesDict, statusText, xhr) {
                    updateExistingChart(coverageScatterMaster, data.coverage.data, 0);
                    updateExistingChart(coverageDetail, newSeriesDict.coverage.data, 0);
                    // updateExistingChart(matchScatterDetail, newSeriesDict.match.data, 0);
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

    });
}