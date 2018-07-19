"use strict";

function move() {
        d3.select(".contain")
            .attr("transform", d3.event.transform);
}

function dataPrep(sankeyData) {
    /*
        Data prep takes the results from the factory promise, and transforms it to the correct format
        for D3-sankey
     */
    let nodes = sankeyData.nodes;
    let sortedLinks = sankeyData.links;
    nodes.links = sankeyData.links.filter(function (link) {
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

    let link = g.append("g")
        .attr("class", "links")
        .attr("fill", "none")
        .attr("stroke", "#000")
        .attr("stroke-opacity", 0.2)
        .selectAll("path");
    let node = g.append("g")
        .attr("class", "nodes")
        .attr("font-family", "sans-serif")
        .attr("font-size", 8)
        .selectAll("g");

    sankey(nodesObj);
    link = link
        .data(nodesObj.links)
        .enter().append("path")
        .attr("d", d3.sankeyLinkHorizontal())
        .attr("class", "link")
        .attr("stroke-width", function (d) {
            return Math.max(1, d.width);
        });

    link.append("title")
        .text(function (d) {
            return `${d.source.name} → ${d.target.name}
${format(d.value)}`;
        });


    node = node
        .data(nodesObj.nodes)
        .enter().append("g").attr("class", "node");

    node.append("rect")
        .attr("x", function (d) {
            return d.x0;
        })
        .attr("y", function (d) {
            return d.y0;
        })
        .attr("height", function (d) {
            return (d.y1 - d.y0);
        })
        .attr("width", function (d) {
            return d.x1 - d.x0;
        })
        .attr("fill", function (d) {
            return color(d.name.replace(/ .*/, ""));
        })
        .attr("stroke", "#000");

    node.append("text")
        .attr("x", function (d) {
            return d.x0 - 6;
        })
        .attr("y", function (d) {
            return (d.y1 + d.y0) / 2;
        })
        .attr("dy", "0.35em")
        .attr("text-anchor", "end")
        .text(function (d) {
            return d.name;
        })
        .filter(function (d) {
            return d.x0 < width / 2;
        })
        .attr("x", function (d) {
            return d.x1 + 6;
        })
        .attr("text-anchor", "start");

    node.append("title")
        .text(function (d) {
            return d.name + "\n" + format(d.value);
        });
}

function update(flowcellId, sankey, checkForData, svg, g, format, color, width) {
    console.log("update");
    $.get("/sankey", {flowcellId}, result => {
        console.log(result);
        let nodes = dataPrep(result.nodes);
        // TODO use is active to check if flowcell is active
        if (result.links.length !== 0 && checkForData) {
            setTimeout(function () {
                d3.select("body").attr("class", "loaded sidebar-collapse")
            }, 1500)
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
    let width = container.getBoundingClientRect().width * 0.957;
    // height of page
    let height = container.getBoundingClientRect().height * 0.6;
    let zoom = d3.zoom()
        .scaleExtent([1, 6]).translateExtent([[0, 0], [width, height * 0.88]])
        .on("zoom", move);
    let svg = d3.select(".sankeyContainer").append("svg").attr("width", width)
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
    console.log(sankey);
    // zoom function, called on svg
    console.log(width, height);

    console.log("Hello there");
    update(flowcellId, sankey, checkForData, svg, g, format, color, width)
}