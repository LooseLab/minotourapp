"use strict";
let updateSankey;
// resize svg and graphic on window resize
// $(window).on("resize", function () {
//     // height of page
//     let hi = $(window).height() * 0.55;
//     // width of page
//     let width = $(window).width() - 60;
//     // flowcell id for data
//     let inputFlowcellId = document.querySelector("#flowcell-id");
//     let flowcellId = inputFlowcellId.value;
//     //update svg width and height
//     d3.select(".svg-sankey").attr("width", width);
//     d3.select(".svg-sankey").attr("height", hi);
//     // draw the sankey diagram
//     drawSankey(flowcellId);
// });

function draw(nodesObj, sankey, g, format, color, width) {
    // Draw the diagram
    let node;
    let text;
    let link;
    let gradient;
    // use d3-sankey to sankeyify the data, providing path coordinates
    sankey(nodesObj);

    // draw the links
    link = g.append("g").attr("class", "links")
        .attr("fill", "none")
        .attr("stroke-opacity", 0.5)
        .selectAll("g")
        .data(nodesObj.links)
        .enter().append("g")
        .style("mix-blend-mode", "multiply");
    // append the linear gradients for the colour
    gradient = link.append("linearGradient")
        .attr("id", function (d, i) {
            return i;
        })
        .attr("gradientUnits", "userSpaceOnUse")
        .attr("x1", d => d.source.x1)
        .attr("x2", d => d.target.x0);
    // set the start colour for the gradient (0%)
    gradient.append("stop")
        .attr("offset", "0%")
        .attr("stop-color", d => color(d.source.name.replace(/ .*/, "")));
    // set the stop colour for the gradient (100%)
    gradient.append("stop")
        .attr("offset", "100%")
        .attr("stop-color", d => color(d.target.name.replace(/ .*/, "")));

    // append the path for the link to rhe SVG
    link.append("path")
        .attr("d", d3.sankeyLinkHorizontal())
        .attr("stroke", function (d, i) {
            return "url(#" + i + ")"
        })
        .attr("stroke-width", d => Math.max(0.5, d.width));
    // append the title
    link.append("title")
        .text(d => `${d.source.name} → ${d.target.name}\n${format(d.value)}`);

    // Append the nodes to the svg
    node = g.append("g").attr("class", "nodes")
        .selectAll("rect")
        .data(nodesObj.nodes)
        .enter().append("rect")
        .attr("x", d => d.x0)
        .attr("y", d => d.y0)
        .attr("height", d => d.y1 - d.y0)
        .attr("width", d => d.x1 - d.x0)
        .attr("fill", function (d) {
            return color(d.name.replace(/ .*/, ""));
        })
        .append("title")
        .text(d => `${d.name}\n${format(d.value)}`);

    // append the text labels to the svg
    text = g.append("g").attr("class", "text")
        .style("font", "10px sans-serif")
        .selectAll("text")
        .data(nodesObj.nodes)
        .enter().append("text")
        .attr("x", d => d.x0 < width / 2 ? d.x1 + 6 : d.x0 - 6)
        .attr("y", d => (d.y1 + d.y0) / 2)
        .attr("dy", "0.35em")
        .attr("text-anchor", d => d.x0 < width / 2 ? "start" : "end")
        .text(d => d.name);

}
// update or draw the existing svg using the AJAX results from the server
function update(flowcellId, sankey, checkForData, svg, g, format, color, width, selectedBarcode) {
    // TODO species limit one day
    $.get("/sankey", {flowcellId , "barcode":selectedBarcode}, result => {
        let nodes;
        // if theres no data from the server
        if (result === undefined) {
            return;
        }

        nodes = result.sankey;
        // If there is data and checkForData is true, clear the css elements and lloading sign so we can see the graphics
        if (result.sankey.nodes.length !== 0 && checkForData) {
            d3.select("#loading-sign").transition().duration(4500).style("opacity", 0);

            setTimeout(function () {
                $("body").addClass("loaded sidebar-collapse");
                d3.select("#loading-sign").style("display", "none");
                d3.select(".vis-container").style("display", "contents");
            }, 5000);
        }
        // TODO update in place
        svg.select(".contain").selectAll("*").remove();
        draw(nodes, sankey, g, format, color, width);
    });
}

function topLevelSankeyDrawer(flowcellID){

    drawSankey(flowcellID);

}

// top level function
function drawSankey(flowcellId) {
    var selectedBarcode = get_selected_barcode();
    // Check the tab value
    let flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');

    if(flowcell_selected_tab_input.value !== "Metagenomics"){
        console.log("clearing snakey");
        clearInterval(updateSankey);
        return;
    }
    //set svg width and height
    // height of page
    let hi = $(window).height() * 0.6;
    // width of page
    hi = d3.max([hi, 480]);
    let width = $(window).width() - 60;
    let svg;
    let g;
    let formatNumber = d3.format(",.0f"),
        format = function (d) {
            return formatNumber(d) + " Reads";
        },
        // create colour scheme
        color = d3.scaleOrdinal(d3.schemeCategory10);
    let checkForData = true;
    let sankey = d3.sankey()
        .nodeWidth(20)
        .nodePadding(3)
        .size([width*0.95, hi * 0.95]).nodeId(function id(d) {
            return d.name;
        })
    ;
    // panning and zoom function
    function move() {
        d3.select(".contain")
            .attr("transform", d3.event.transform);
    }

    let zoom = d3.zoom()
        .scaleExtent([1, 6]).translateExtent([[0, 0], [width, hi]])
        .on("zoom", move);
    if ($(".svg-sankey").length !== 0) {
        svg = d3.select(".svg-sankey");
        g = d3.select(".contain");
    } else {
        console.log(hi);
        svg = d3.select(".svg-container").append("svg").attr("width", width).attr("class", "svg-sankey")
            .attr("height", hi).call(zoom);
        g = svg.append("g").attr("class", "contain");
    }
    update(flowcellId, sankey, checkForData, svg, g, format, color, width, selectedBarcode);
    let timmy = setTimeout(drawSankey, 60000, flowcellId, selectedBarcode);
}
