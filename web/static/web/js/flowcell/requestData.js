function requestData(flowcell_id) {

    var selected_barcode = get_selected_barcode();

    var flowcell_id = flowcell_id;

    var url_run = '/api/v1/flowcells/' + flowcell_id + '/';

    $.get(url_run, (function (data) {

        var barcodes = new Set();

        for (var i = 0; i < data.barcodes.length; i++) {
            barcodes.add(data.barcodes[i].name);
        }

        var flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');

        if (flowcell_selected_tab_input.value == 'Summary') {

            this.requestRunDetails(flowcell_id);
            requestMinknowMessages(flowcell_id, data);

        } else if (flowcell_selected_tab_input.value == 'Tasks') {

            this.flowcellTaskHistoryTable(flowcell_id);

        } else if (flowcell_selected_tab_input.value == 'Basecalled Data') {

            if (selected_barcode == '') {
                set_selected_barcode('All reads');
            }

            this.barcodes = Array.from(barcodes).sort();
            this.updateBarcodeNavTab();
            this.requestSummaryData(flowcell_id);
            this.requestHistogramData(flowcell_id);
            this.requestStatistics(flowcell_id);
            this.requestChannelSummaryData(flowcell_id);
            //this.requestReference(flowcell_id);

        } else if (flowcell_selected_tab_input.value == 'Live Event Data') {

            this.requestRunDetails(flowcell_id);
            this.requestLiveRunStats(flowcell_id);

        } else if (flowcell_selected_tab_input.value == 'Runs') {

        } else if (flowcell_selected_tab_input.value == 'Metagenomics') {
            // The intervals for updating the charts are found in the individual files in the vis-d3 directory
            // SO you are on the Metagenomics tab
            if (selected_barcode == '') {
                set_selected_barcode('All reads');
            }

            this.barcodes = Array.from(barcodes).sort();
            // add the barcodes for this flowcell
            this.addBarcodeTabs(flowcell_id);
            // draw the sankey
            this.topLevelSankeyDrawer(flowcell_id);
            // update the metadata header
            this.topMetaHeader(flowcell_id);
            // Draw the donut chart
            this.topLevelDrawDonut(flowcell_id);
            // update the total Reads Table
            this.topGetTotalReadsTable(flowcell_id);
            // Draw the donut rank table
            this.topGetDonutRankTable(flowcellId);
            // Draw the alert mapping targets table;
            this.update_mapping_table(flowcellId);

        } else if (flowcell_selected_tab_input.value == 'Sequence Mapping') {

            this.requestPafData(flowcell_id);

        } else if (flowcell_selected_tab_input.value == 'Assembly') {

            this.requestGfaData(flowcell_id);

        }
    }).bind(this));
    setTimeout(requestData.bind(this, flowcell_id), 60000);
};
