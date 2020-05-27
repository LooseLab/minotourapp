function requestData(flowcell_id) {
    console.log("Requesting data");
    var selected_barcode = getSelectedBarcode("Metagenomics");
    console.log(selected_barcode)
    var flowcell_id = flowcell_id;

    var url_run = '/api/v1/flowcells/' + flowcell_id + '/';

    $.get(url_run, (function (result) {
        // Get data about teflowcell including barcode and number of runs
        let data = result.data;

        var barcodes = new Set();
        // For each of the barcodes, add it to the set to be used to draw and select the barcode tabs later
        for (var i = 0; i < data.barcodes.length; i++) {
            if (!['No barcode','S','U'].includes(data.barcodes[i].name)) {
                barcodes.add(data.barcodes[i].name);
            }
        }

        var flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');

        if (flowcell_selected_tab_input.value == 'nav-summary-data') {
            // if we are on the summary page - cal l requestRunDetail from the monitor app scope
            // appends the Run summary table HTML onto the summary page
            this.requestRunDetails(flowcell_id);
            //// append a messages html tableonto the messages div of the summaries tab
            requestMinknowMessages(flowcell_id, data);

        } else if (flowcell_selected_tab_input.value == 'nav-tasks') {
            // If we are on the tasks tab
            // Load the tasks form for selecting a task and reference and submitting a ne wtask
            loadTasksForm();
            // Load the flowcell tasks history table

            this.flowcellTaskHistoryTable(flowcell_id);

        } else if (flowcell_selected_tab_input.value == 'nav-reads') {

            requestReadData(flowcell_id);

        } else if (flowcell_selected_tab_input.value == 'nav-live-event-data') {

            // this.requestRunDetails(flowcell_id);
            this.requestLiveRunStats(flowcell_id);

        } else if (flowcell_selected_tab_input.value == 'nav-metagenomics') {
            // The intervals for updating the charts are found in the individual files in the vis-d3 directory
            // SO you are on the Metagenomics tab, congratulations!
            let addBarcodes = this.addMetaBarcodeTabs.bind(this);
            if (selected_barcode == null) {
                setSelectedBarcode('All reads', "Metagenomics");
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
            // this.requestAdvancedPafData(flowcell_id);
            this.drawReadUntilCharts();

        } else if (flowcell_selected_tab_input.value == 'nav-sequence-assembly') {

            this.requestGfaData(flowcell_id);

        } else if (flowcell_selected_tab_input.value == "nav-notifications") {
            notificationsController.getReferencesForDataList();
        }

    }).bind(this));
};
