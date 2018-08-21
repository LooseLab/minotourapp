"use strict";

function move() {
    d3.select(".badCopNoDonut").attr("transform", d3.event.transform);
}

function drawPie(countedData, dataLength, fillArray, pie, arc, svg) {
    // draw the dount chart
    let color = d3.scaleOrdinal(d3.schemeCategory10);
    let slice = svg.select(".slices").selectAll("path.slice")
        .data(pie(countedData));
    //enter your data, returned from pie(countedData, insert a path, set d to data provided by arc function)
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

    slice.exit().remove();
}

function drawDonut(flowCellId) {
    let running = true;
    let taxas = ["species", "genus", "family", "order", "classy", "phylum", "superkingdom"];
    // the taxa titles we wish to display under the slider
    let DisplayTaxas = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom"];
    let container = d3.select("body").node();
    let width = ((container.getBoundingClientRect().width - 30) * 0.25) - 20;
    let height = $(window).height() * 0.35,
        radius = Math.min(width, height / 1.7);
    let zoom = d3.zoom()
        .scaleExtent([1, 10]).translateExtent([[0, 0], [width, height]])
        .on("zoom", move);
    let number;
    let range = $('.input-range'),
        value = $('.taxa-level');
    let dHeight = 35;
    let svg;
    let g;
    // pie is a d3 function that transforms the value for each of the .values in the array of objects, into
    let pie = d3.pie().padAngle(0.01).sort(null)
        .value(function (d) {
            return d.value;
        });
    // arc is a d3 function that generates arc paths from the data provided from d3.pie
    let arc = d3.arc()
        .innerRadius(radius * 0.8)
        .outerRadius(radius * 0.5)
        .padAngle(0.02);

    const fillArray = ["rgb(255, 0, 0)", "rgb(0, 0, 255)", "rgb(0, 255, 0)", "rgb(255, 0, 255)", "rgb(255, 153, 51)", "rgb(255, 252, 8)",
        "rgb(8, 249, 255)", "rgb(202, 8 ,255)", "rgb(35, 144, 35)", "rgb(202, 119, 8)", "rgb(173, 176, 150)", "rgb(0, 0, 0)",
        "rgb(141, 202, 8)", "rgb(8, 204, 200)", "rgb(117, 114, 45)", "rgb(255, 0, 99)", "rgb(114, 132, 114)", "rgb(204, 255, 102)",
        "rgb(204, 153, 255)", "rgb(255, 153, 0)"];

    if ($(".donut-svg").length) {
        svg = d3.select(".donut-svg");
        // d3.select(".slices")/**/.selectAll("*").remove();
    } else {
        svg = d3.select(".donutContainer").append("svg")
            .attr("class", "donut-svg")
            .attr("width", width)
            .style("height", dHeight+"vh")
            .attr("margin", "auto")
            .style("margin-left", "-20px")
            .attr("display", "block")
            .call(zoom);
        g = svg.append("g").attr("class", "badCopNoDonut");
        // append g elements to the svg, attr transform translate centers them in the svg, otherwise drawn offscreen

        g.append("g")
            .attr("class", "slices")
            .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");
    }

    value.html(DisplayTaxas[0]);
    $.get("/donut", {flowcellId: flowCellId, visType: "donut"}, result => {
        console.log(result)
        let dataToDraw = result.result;
        // if there is no data return and try again when interval is up on $interval
        if (result.result.length === 0) {
            return;
        }
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
            // remove the current donut chart
            // d3.select(".slices").selectAll("*").remove();
            // get the right taxa for the key to the results object
            let current_selected_taxa = taxas[number];
            // get the results data array for the currently selected taxa clade
            let sortedData = dataToDraw[number][current_selected_taxa];
            // datalength - how many members ar ein this clade (1-20)
            dataLength = sortedData.length;
            //draw a new donut
            drawPie(sortedData, dataLength, fillArray, pie, arc, svg);
        });
        drawPie(dataToDraw[0]["species"], dataLength, fillArray, pie, arc, svg);
    });
}
