function requestData(flowcell_id) {

    var selected_barcode = get_selected_barcode();

    var flowcellId = flowcell_id;

    var url_run = '/api/v1/flowcells/' + flowcellId + '/';

    $.get(url_run, (function (data) {

        var barcodes = new Set();

        for (var i = 0; i < data.barcodes.length; i++) {
            barcodes.add(data.barcodes[i].name);
        }

        var flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');

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
            // The intervals for updating the charts are found in the individual files in the vis-d3 directory
            // SO you are on the Metagenomics tab
            if (selected_barcode == '') {
                set_selected_barcode('All reads');
            }

            this.barcodes = Array.from(barcodes).sort();
            // add the barcodes for this flowcell
            this.addBarcodeTabs(flowcellId);
            // draw the sankey
            this.topLevelSankeyDrawer(flowcellId, selected_barcode);
            // update the metadata header
            this.topMetaHeader(flowcellId);
            // Draw the donut chart
            this.topLevelDrawDonut(flowcellId, selected_barcode);
            // update the total Reads Table
            this.topGetTotalReadsTable(flowcellId, selected_barcode);
            // Draw rhe donut rank table
            this.topGetDonutRankTable(flowcellId, selected_barcode);

        } else if (flowcell_selected_tab_input.value == 'Sequence Mapping') {

            this.requestPafData(flowcellId);

        } else if (flowcell_selected_tab_input.value == 'Assembly') {

            this.requestGfaData(flowcellId);

        }
    }).bind(this));
    // todo is this fixed?
    // setTimeout(requestData.bind(this, flowcell_id), 60000);
};
