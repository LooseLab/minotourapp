let comparison_tab = {};
function requestData(flowcell_id) {
    let updateSankey;
    let updateMeatheader;
    let updateDonut;
    let updateTotalTable;
    let updateDonutTable;
    let check;
    function clearIntervals(){
        clearInterval(updateSankey);
        clearInterval(updateMeatheader);
        clearInterval(updateDonut);
        clearInterval(updateTotalTable);
        clearInterval(updateDonutTable);
        clearInterval(check);
        console.log("clearing intervals");
    }
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
        if (comparison_tab.tab_value !== flowcell_selected_tab_input.value) {
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
            comparison_tab.tab_value = flowcell_selected_tab_input.value;
            clearIntervals();
            checkRunning(flowcellId);
            this.drawSankey(flowcellId);
            updateSankey = setInterval(this.drawSankey, 60000, flowcellId);
            this.metaHeader(flowcellId);
            updateMeatheader = setInterval(this.metaHeader, 60000, flowcellId);
            this.drawDonut(flowcellId);
            updateDonut = setInterval(this.drawDonut, 60000, flowcellId);
            this.getTotalReadsTable(flowcellId);
            updateTotalTable = setInterval(this.getTotalReadsTable, 60000, flowcellId);
            this.getDonutRankTable(flowcellId);
            updateDonutTable = setInterval(this.getDonutRankTable, 60000, flowcellId);
            checkRunning(flowcellId);
            check = setInterval(checkRunning, 30000, flowcellId);
            $(window).on("unload", function () {
                clearInterval(updateSankey);
                clearInterval(updateMeatheader);
                clearInterval(updateDonut);
                clearInterval(updateTotalTable);
                clearInterval(updateDonutTable);
                clearInterval(check);
            });

        } else if (flowcell_selected_tab_input.value == 'Sequence Mapping') {

            this.requestPafData(flowcellId);

        } else if (flowcell_selected_tab_input.value == 'Assembly') {

            this.requestGfaData(flowcellId);

        }
    }).bind(this));

    // setTimeout(requestData.bind(this, flowcell_id), 60000);
};
