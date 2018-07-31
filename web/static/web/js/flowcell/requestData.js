function requestData(flowcell_id) {

    console.log('Running requestData for flowcell ' + flowcell_id);

    var selected_barcode = get_selected_barcode();

    var flowcellId = flowcell_id;

    var url_run = '/api/v1/flowcells/' + flowcellId + '/';

    $.get(url_run, (function (data) {

        var barcodes = new Set();

        for (var i = 0; i < data.barcodes.length; i++) {
            barcodes.add(data.barcodes[i].name)
        }

        var flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');

        this.rundata = data;  // TODO who is using this variable?

        if (flowcell_selected_tab_input.value == 'Summary') {

            this.requestRunDetails(flowcellId);
            requestMinknowMessages(flowcellId, data);

        } else if (flowcell_selected_tab_input.value == 'Tasks') {

            this.requestTasks(flowcellId);
            this.liveUpdateTasks(flowcellId);

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
            this.drawSankey(flowcellId);
            this.metaHeader(flowcellId);
            this.drawDonut(flowcellId);
            this.getTotalReadsTable(flowcellId);
            this.getDonutRankTable(flowcellId);

        } else if (flowcell_selected_tab_input.value == 'Sequence Mapping') {

            this.requestPafData(flowcellId);

        } else if (flowcell_selected_tab_input.value == 'Assembly') {

            this.requestGfaData(flowcellId);

        }
    }).bind(this));

    // setTimeout(requestData.bind(this, flowcell_id), 60000);
};
