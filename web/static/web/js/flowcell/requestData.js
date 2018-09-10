let comparison_tab = "";
// Comparison tab stores the
function requestData(flowcell_id) {
    let updateSankey;
    let updateMeatheader;
    let updateDonut;
    let updateTotalTable;
    let updateDonutTable;
    let check;

    // Clear interavls clears all the interbal calls that we are using to update the page
    function clearIntervals(){
        clearInterval(updateSankey);
        clearInterval(updateMeatheader);
        clearInterval(updateDonut);
        clearInterval(updateTotalTable);
        clearInterval(updateDonutTable);
        clearInterval(check);
        console.log("clearing intervals");
    }
    // Check if the Flowcell is still running, if it isn't cnacel the intervals
    function checkRunning(flowcellID){
        let url = '/api/v1/tasks/';
        let params = {"search_criteria": "flowcell", "search_value": flowcellID};
        $.get(url, params, function (result) {
            let data = result.data;
            data.sort((a, b) => b.id - a.id);
            console.log(data);
            if(data[0].complete){
                console.log("complete");
                clearIntervals();
            }
        });
        console.log("check running");
    }

    console.log('Running requestData for flowcell ' + flowcell_id);



    var selected_barcode = get_selected_barcode();

    var flowcellId = flowcell_id;

    var url_run = '/api/v1/flowcells/' + flowcellId + '/';

    $.get(url_run, (function (data) {

        var barcodes = new Set();

        for (var i = 0; i < data.barcodes.length; i++) {
            barcodes.add(data.barcodes[i].name);
        }

        var flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');
        // If you've changed tabs, clear the intervals that are updating the metagenomics visualisations
        if (comparison_tab !== flowcell_selected_tab_input.value) {
            console.log("changed tabs, clearing intervals");
            clearIntervals();
        }

        this.rundata = data;  // TODO who is using this variable?
        // TODO I don't know man

        if (flowcell_selected_tab_input.value == 'Summary') {

            this.requestRunDetails(flowcellId);
            requestMinknowMessages(flowcellId, data);

        } else if (flowcell_selected_tab_input.value == 'Tasks') {
            this.flowcellTaskHistoryTable(flowcellId);

        } else if (flowcell_selected_tab_input.value == 'Basecalled Data') {
            if (selected_barcode == '') {
                set_selected_barcode('All reads');
            }

            this.barcodes = Array.from(barcodes).sort();
            this.updateBarcodeNavTab();
            this.requestSummaryData(flowcellId);
            this.requestHistogramData(flowcellId);
            this.requestStatistics(flowcellId);
            this.requestChannelSummaryData(flowcellId);
            //this.requestReference(flowcellId);

        } else if (flowcell_selected_tab_input.value == 'Live Event Data') {

            this.requestRunDetails(flowcellId);
            this.requestLiveRunStats(flowcellId);

        } else if (flowcell_selected_tab_input.value == 'Runs') {

        } else if (flowcell_selected_tab_input.value == 'Metagenomics') {
            // TODO rewrite into individual funcitons, find a way to note href change
            // SO you are on the Metagenomics tab
            comparison_tab = flowcell_selected_tab_input.value;
            // Clear any previous intervals if they've somehow survived
            clearIntervals();
            // Cheeck to see if the flowcell Job is still running
            checkRunning(flowcellId);
            // Draw the sankey diagram
            this.drawSankey(flowcellId);
            // Setup the interval request fto update the sankey
            updateSankey = setInterval(this.drawSankey, 60000, flowcellId);
            // update the metadata header
            this.metaHeader(flowcellId);
            // Setup metaheader interval
            updateMeatheader = setInterval(this.metaHeader, 60000, flowcellId);
            // Draw the donut chart
            this.drawDonut(flowcellId);
            // Update the Donut Chart by creating an interval for every 60 seconds
            updateDonut = setInterval(this.drawDonut, 60000, flowcellId);
            // update the total Reads Table
            this.getTotalReadsTable(flowcellId);
            // Create a 60 secons unterval to update the table
            updateTotalTable = setInterval(this.getTotalReadsTable, 60000, flowcellId);
            // Draw rhe donut rank table
            this.getDonutRankTable(flowcellId);
            // Create 60 second update for donut table
            updateDonutTable = setInterval(this.getDonutRankTable, 60000, flowcellId);
            // Check to see if the analaysis is still running
            checkRunning(flowcellId);
            // 30 second interval for flowcell is still running
            check = setInterval(checkRunning, 30000, flowcellId);
            // If the page changes, clear the intervals
            $(window).on("unload", clearIntervals);

        } else if (flowcell_selected_tab_input.value == 'Sequence Mapping') {

            this.requestPafData(flowcellId);

        } else if (flowcell_selected_tab_input.value == 'Assembly') {

            this.requestGfaData(flowcellId);

        }
    }).bind(this));

    // setTimeout(requestData.bind(this, flowcell_id), 60000);
};
