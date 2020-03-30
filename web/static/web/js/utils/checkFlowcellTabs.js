// function checkFlowcellTabs(flowcell_id) {
//
//     var url = "/api/v1/flowcells/" + flowcell_id + "/tabs/";
//
//     $.getJSON(url, function (data) {
//
//         var items = [];
//
//         $.each(data, function (key, value) {
//             items.push(value);
//         });
//
//         items.sort(function (a, b) {
//             return a.position - b.position;
//         });
//
//         var flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');
//     });
//
//     setTimeout(checkFlowcellTabs, 60000, flowcell_id);
// }

