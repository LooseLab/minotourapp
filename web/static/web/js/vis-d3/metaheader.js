"use strict";
let firsty = true;
function metaHeader(flowcellId){
    let table, head, row;
    if(firsty){
        table = d3.select(".meta_taberu").append("table").attr("class", "table table-hover").attr("table-layout", "fixed");
        head = table.append("thead");
        row = head.append("tr");
    } else{
        row = d3.select(".meta_taberu").select("table").select("tr");
    }
    $.get("/metaview", {flowcellId}, result => {
        if(result[1].value === 0){
            return;
        }
        row.selectAll("th").data(result).enter().append("th");

        row.selectAll("th").data(result).html(function(d){
            return d.key + d.value;
        });
        }
    ).fail(function(){
        return;
    });
    firsty = false;
}