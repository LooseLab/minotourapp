"use strict";
function move() {
        d3.select(".badCopNoDonut").attr("transform", d3.event.transform);
}
function setup() {
    let running = true;
    let DisplayTaxas = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subspecies", "Strain"];
    let taxas = ["superkingdom", "phylum", "classy", "order", "family", "genus", "species", "subspecies", "strain"];
    let container = d3.select(".donutContainer").node();
    let width = container.getBoundingClientRect().width;
    let height = container.getBoundingClientRect().height * 1.5,
        radius = Math.min(width, height / 1.7);
    let zoom = d3.zoom()
        .scaleExtent([1, 10]).translateExtent([[0, 0], [width, height]])
        .on("zoom", move);
    let number;
    let range = $('.input-range'),
        value = $('.taxa-level');
// reference the d3 selection of the svg element
    let svg = d3.select(".donut")
        .attr("width", width)
        .attr("height", height)
        .attr("margin", "auto")
        .style("margin-left", "-20px")
        .attr("display", "block")
        .call(zoom);
// set the duration of the transition
    let duration = 3000;
// append g elements to the svg, attr transform translate centers them in the svg, otherwise drawn offscreen
    let g = svg.append("g").attr("class", "badCopNoDonut");
    // pie is a d3 function that transforms the value for each of the .values in the array of objects, into
    let pie = d3.pie().padAngle(0.01).sort(null)
        .value(function (d) {
            return d.value;
        });
    // arc is a d3 function that generates arc paths from the data provided from d3.pie
    let arc = d3.arc()
        .innerRadius(radius * 0.8)
        .outerRadius(radius * 0.5)
        .padAngle(0.05);
    g.append("g")
        .attr("class", "slices")
        .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");
    value.html(DisplayTaxas[0]);
    $.get("/donut", {meta_id: 5}, result => {
        console.log(result);
    });

}
