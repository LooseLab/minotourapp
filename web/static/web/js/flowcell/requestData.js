function requestData(flowcell_id) {

    var selected_barcode = get_selected_barcode();

    var flowcell_id = flowcell_id;

    var url_run = '/api/v1/flowcells/' + flowcell_id + '/';

    $.get(url_run, (function (result) {

        let data = result.data;
        // let meta_barcodes= result.meta_barcodes;
        var barcodes = new Set();

        for (var i = 0; i < data.barcodes.length; i++) {
            if (data.barcodes[i].name != 'No barcode') {
                barcodes.add(data.barcodes[i].name);
            }
        }

        var flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');

        if (flowcell_selected_tab_input.value == 'nav-summary-data') {

            this.requestRunDetails(flowcell_id);
            console.log(data);
            requestMinknowMessages(flowcell_id, data);

        } else if (flowcell_selected_tab_input.value == 'nav-tasks') {

            this.flowcellTaskHistoryTable(flowcell_id);

        } else if (flowcell_selected_tab_input.value == 'nav-basecalled-data') {

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

        } else if (flowcell_selected_tab_input.value == 'nav-reads') {

            requestReadData(flowcell_id);

        } else if (flowcell_selected_tab_input.value == 'nav-live-event-data') {

            // this.requestRunDetails(flowcell_id);
            this.requestLiveRunStats(flowcell_id);

        } else if (flowcell_selected_tab_input.value == 'nav-metagenomics') {
            // The intervals for updating the charts are found in the individual files in the vis-d3 directory
            // SO you are on the Metagenomics tab, congratulations!
            let addBarcodes = this.addMetaBarcodeTabs.bind(this);
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
            // draw the simplified results table
            this.draw_simple_table(flowcell_id);

            $.get(url, {}, function (result) {
                this.barcodes = result.data.sort();

                this.tabs = result.tabs;
                // draw the sankey
                addBarcodes(flowcell_id, this.barcodes, this.tabs);
                // update the metadata header

            });


        } else if (flowcell_selected_tab_input.value == 'nav-sequence-mapping') {
            this.requestPafData(flowcell_id);

        } else if (flowcell_selected_tab_input.value == 'nav-advanced-sequence-mapping') {
            this.requestAdvancedPafData(flowcell_id);
            this.drawReadUntilCharts();

        } else if (flowcell_selected_tab_input.value == 'nav-sequence-assembly') {

            this.requestGfaData(flowcell_id);

        }
    }).bind(this));
};
