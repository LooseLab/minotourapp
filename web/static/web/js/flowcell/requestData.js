function requestData(flowcell_id) {

    var selected_barcode = get_selected_barcode();

    var flowcell_id = flowcell_id;

    var url_run = '/api/v1/flowcells/' + flowcell_id + '/';
    $.get(url_run, (function (result) {

        let data = result.data;
        // let meta_barcodes= result.meta_barcodes;
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
            let bants = this.addMetaBarcodeTabs.bind(this);
            if (selected_barcode == '') {
                set_selected_barcode('All reads');
            }
            let url = "/api/v1/flowcells/" + flowcell_id + "/metagenomic_barcodes";
            this.drawSankey(flowcell_id);
            this.metaHeader(flowcell_id);
            // Draw the donut rank table
            this.drawDonutRankTable(flowcell_id);
            // Draw the donut chart
            this.drawDonut(flowcell_id);
            // update the total Reads Table
            this.getTotalReadsTable(flowcell_id);
            // Draw the alert mapping targets table;
            this.update_mapping_table(flowcell_id);

            $.get(url, {}, function (result) {
                this.barcodes = result.data.sort();
                // draw the sankey
                bants(flowcell_id, this.barcodes);
                // update the metadata header

            });

        } else if (flowcell_selected_tab_input.value == 'Sequence Mapping') {

            this.requestPafData(flowcell_id);

        } else if (flowcell_selected_tab_input.value == 'Assembly') {

            this.requestGfaData(flowcell_id);

        }
    }).bind(this));
    // setTimeout(requestData.bind(this, flowcell_id), 60000);
};
