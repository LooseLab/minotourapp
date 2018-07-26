"use strict";


function dataPrep(sankeyData) {
    /*
        Data prep takes the results from the factory promise, and transforms it to the correct format
        for D3-sankey
     */

    let nodes = sankeyData.queryset;
    let sortedLinks = sankeyData.values;
    nodes.links = nodes.links.filter(function (link) {
        return link.target !== link.source;
    });
    let nodeNames = [];
    let topFiftySpecies = sortedLinks.filter(function (link) {
        return link.target_tax_level === "species";
    }).splice(0, 50);
    let topFiftyLink = [];

    for (let i = 0; i < topFiftySpecies.length; i++) {
        let specimin = topFiftySpecies[i];
        topFiftyLink.push(specimin);
        nodeNames.push({"name": specimin.source}, {"name": specimin.target});
        let notTopOfTree = true;
        let first = true;
        while (notTopOfTree) {
            let source;
            let link;
            if (first) {
                source = specimin.source;
                first = false;
            } else {
                let sourcey = topFiftyLink[topFiftyLink.length - 1];
                source = sourcey.source;
            }
            if (nodes.links.findIndex(x => x.target === source) ===
                -1
            ) {
                notTopOfTree = false;
            }
            if (notTopOfTree) {
                link = nodes.links.filter(function (l) {
                    return l.target === source
                });
                topFiftyLink.push(link[0]);
                nodeNames.push({"name": link[0].target}, {"name": link[0].source})
            }
        }
    }
    nodeNames = nodeNames.filter((node, index, self) => self.findIndex(t => t.name === node.name
        ) ===
        index
    )
    ;
    topFiftyLink = topFiftyLink.filter((link, index, self) => self.findIndex(t => t.source === link.source && t.target === link.target
        ) ===
        index
    )
    ;
    nodes.links = topFiftyLink;
    nodes.nodes = nodeNames;
    return nodes;
}

function draw(nodesObj, sankey, g, format, color, width) {
    sankey(nodesObj);

    const link = g.append("g").attr("class", "links")
        .attr("fill", "none")
        .attr("stroke-opacity", 0.5)
        .selectAll("g")
        .data(nodesObj.links)
        .enter().append("g")
        .style("mix-blend-mode", "multiply");

      const gradient = link.append("linearGradient")
      .attr("id", function(d, i){
          return i;
      })
      .attr("gradientUnits", "userSpaceOnUse")
      .attr("x1", d => d.source.x1)
      .attr("x2", d => d.target.x0);

      gradient.append("stop")
          .attr("offset", "0%")
          .attr("stop-color", d => color(d.source.name.replace(/ .*/, "")));

      gradient.append("stop")
          .attr("offset", "100%")
          .attr("stop-color", d => color(d.target.name.replace(/ .*/, "")));

    link.append("path")
        .attr("d", d3.sankeyLinkHorizontal())
        .attr("stroke", function(d, i){
           return "url(#" + i+ ")"
         })
        .attr("stroke-width", d => Math.max(1, d.width));

    link.append("title")
      .text(d => `${d.source.name} → ${d.target.name}\n${format(d.value)}`);

    let node = g.append("g").attr("class", "nodes")
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

    let text = g.append("g").attr("class", "text")
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

function update(flowcellId, sankey, checkForData, svg, g, format, color, width) {
    console.log("update");
    $.get("/sankey", {flowcellId}, result => {
        console.log(result);
        let nodes = dataPrep(result);
        console.log(nodes);
        // TODO use is active to check if flowcell is active
        if (result.queryset.nodes.length !== 0 && checkForData) {
            setTimeout(function () {
                d3.select("body").attr("class", "loaded sidebar-collapse")
            }, 1500);
            checkForData = false;
        }
        svg.select(".contain").selectAll("*").remove();
        draw(nodes, sankey, g, format, color, width);
    });
}

function testHello(flowcellId) {
    let running = true;
    //set svg width and height
    let container = d3.select("body").node();
    // width of page
    let width = container.getBoundingClientRect().width * 0.933;
    // height of page
    let height = container.getBoundingClientRect().height * 0.6;

    function move() {
        d3.select(".contain")
            .attr("transform", d3.event.transform);
    }

    let zoom = d3.zoom()
        .scaleExtent([1, 6]).translateExtent([[0, 0], [width, height * 0.9]])
        .on("zoom", move);
    let svg = d3.select(".svg-container").append("svg").attr("width", width)
        .attr("height", height * 0.9).call(zoom);

    let formatNumber = d3.format(",.0f"),
        format = function (d) {
            return formatNumber(d) + " Reads";
        },
        // create colour scheme
        color = d3.scaleOrdinal(d3.schemeCategory10);
    let g = svg.append("g").attr("class", "contain");
    let checkForData = true;
    let sankey = d3.sankey()
        .nodeWidth(20)
        .nodePadding(5)
        .size([width, height * 0.88]).nodeId(function id(d) {
            return d.name;
        })
    ;
    update(flowcellId, sankey, checkForData, svg, g, format, color, width)
}