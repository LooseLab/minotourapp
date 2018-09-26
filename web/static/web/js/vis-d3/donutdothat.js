"use strict";
// Redraw the SVGs on window resize
let updateDonut;

// $(window).on("resize", function(){
//     let width = ($(window).width() * 0.25) - 50;
//     let height = $(window).height() * 0.35;
//     // Update the svg width and height
//     d3.select(".donut-svg").attr("width", width);
//     d3.select(".donut-svg").attr("height", height);
//     // Get the flowcell id
//     let inputFlowcellId = document.querySelector("#flowcell-id");
//     let flowcellId = inputFlowcellId.value;
//     // Redraw the donut
//     drawDonut(flowcellId);
// });

// The panning and zooming function, called when you apply a call of zoom to the svg on initialisation
function move() {
    d3.select(".badCopNoDonut").attr("transform", d3.event.transform);
}

// draw the actual donut chart slices
function drawPie(countedData, pie, arc, svg) {
    // draw the dount chart
    // select an actual good colour scheme
    let color = d3.scaleOrdinal(d3.schemeCategory10);
    // Select all the slices and bind the provided data after it's been transformed by d3.pie()
    let slice = svg.select(".slices").selectAll("path.slice")
        .data(pie(countedData));
    //enter your data, returned from pie(countedData, insert a path, set d to data provided by arc function)
    // Update existing donut slices
    slice.attr("class", "slice")
        .style("fill", function (d) {
            return color(d.data.label);
        })
        .attr("stroke", "black")
        .attr("stroke-width", 0.2)
        .attr("d", arc).select("title")
        .text(function (d) {
            return d.data.label + "\n" + d.value + " reads";
        });
    // add the new donut slices, held in d3s enter selection
    slice.enter()
        .insert("path")
        .attr("class", "slice")
        .style("fill", function (d, i) {
            return color(d.data.label);
        })
        .attr("stroke", "black")
        .attr("stroke-width", 0.2)
        .attr("d", arc).append("title")
        .text(function (d) {
            return d.data.label + "\n" + d.value + " reads";
        });
    // remove any slices held in the exit selection, that used to have DOM elements but now have no data for them
    slice.exit().remove();
}
// TODO this could be more effeicent as it is sometimes called unnecessarily;
function topLevelDrawDonut(flowCellId, selectedBarcode){
    drawDonut(flowCellId, selectedBarcode);
    updateDonut = setInterval(drawDonut, 60000, flowCellId);
}

function drawDonut(flowCellId, selectedBarcode) {
    let flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');
    if(flowcell_selected_tab_input.value !== "Metagenomics"){
        clearInterval(updateDonut);
        console.log("cleared donut interval");
        return;
    }
    // setup the donut chart
    // the taxas in the order we want, to access them from the AJAX get request results
    let taxas = ["species", "genus", "family", "order", "classy", "phylum", "superkingdom"];
    // the taxa titles we wish to display under the slider
    let DisplayTaxas = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom"];
    // Calculate the width
    let width = ($(window).width() *0.25) -50;
    // Calculate the height
    let height = $(window).height() * 0.35,
        // the radius, the smallest of the width and height /2 so it fits in hte svg
        radius = Math.min(width, height)/2;
    // The d3 zoom function
    let zoom = d3.zoom()
        .scaleExtent([1, 10]).translateExtent([[0, 0], [width, height]])
        .on("zoom", move);
    let number;
    // Select the Range of the html slider and the Displayed label underneath the slider
    let range = $('.input-range'),
        value = $('.taxa-level');
    // Declare the svg and group element variables
    let svg;
    let g;
    // pie is a d3 function that transforms the value for each of the .values in the array of objects, into
    let pie = d3.pie().padAngle(0.01).sort(null)
        .value(function (d) {
            return d.value;
        });
    // arc is a d3 function that generates arc paths from the data provided from d3.pie
    let arc = d3.arc()
        // \The inner and outer radius of the actual donut slices
        .innerRadius(radius * 0.8)
        .outerRadius(radius * 0.5)
        .padAngle(0.02);

    // If there is already a donut-svg
    if ($(".donut-svg").length !== 0) {
        // Select the svg
        svg = d3.select(".donut-svg");
        // recenter, as resizing doesn't change the tansformation
        d3.select(".slices").attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");
    } else {
        // If no extant svg, append a new svg and g element
        svg = d3.select(".donutContainer").append("svg")
            .attr("class", "donut-svg")
            .attr("width", width)
            .attr("height", height)
            .attr("margin", "auto")
            .attr("display", "block")
            .call(zoom);
        g = svg.append("g").attr("class", "badCopNoDonut");

        // append g elements to the svg, attr transform translate centers them in the svg, otherwise drawn offscreen
        g.append("g")
            .attr("class", "slices")
            .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");
        // Set the slider level to Display Species
        value.html(DisplayTaxas[0]);
    }
    // Get the data from the server to sisplay
    $.get("/donut", {flowcellId: flowCellId, visType: "donut", barcode: selectedBarcode}, result => {
        // if there is no data return and try again when interval is up on $interval
        if (result === undefined) {
            return;
        }
        console.log(result);
        let dataToDraw = result.result;
        // what the label is displaying, starts on species
        let currentSelectionSlider = value.html();
        // what index is that in the display taxas so we can get the right value from the internal taxas array and results array
        let index = DisplayTaxas.indexOf(currentSelectionSlider);
        // get the right taxa string, so we can use it as a key o the results object
        let currently_selected_taxa = taxas[index];
        // data1 is the data for the currently selected taxa from the results array
        let data1 = dataToDraw[index][currently_selected_taxa];
        // get the number of members in this clade used to determine whether we need a second table
        let dataLength = data1.length;
        // if the range slider is changed call the anonymous function to redraw everything
        range.on("input", function () {
            // the selected number for the slider level (0-6)
            number = this.value;
            // set the html below the slider to the right level
            value.html(DisplayTaxas[number]);
            // get the right taxa for the key to the results object
            let current_selected_taxa = taxas[number];
            // get the results data array for the currently selected taxa clade
            let sortedData = dataToDraw[number][current_selected_taxa];
            // datalength - how many members ar ein this clade (1-20)
            dataLength = sortedData.length;
            //draw a new donut
            drawPie(sortedData, pie, arc, svg);
        });
        // Draw a new pie chart
        drawPie(data1, pie, arc, svg);
    });
}
