function requestData(flowcell_id) {

    var selected_barcode = get_selected_barcode();

    var flowcellId = flowcell_id;

    var url_run = '/api/v1/flowcells/' + flowcellId + '/';

    $.get(url_run, (function (data) {

        if (this.summary !== null) {
            if (this.summary["All reads"]['Template']["max_length"]["all"]["data"][0] >= 1000000 && this.millionaire != true) {
                $('#eastermodal').modal('show');
                this.millionaire = true;
            }
        }

        var barcodes = new Set();

        for (var i = 0; i < data.barcodes.length; i++) {
            barcodes.add(data.barcodes[i].name)
        }

        var flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');

        if (flowcell_selected_tab_input.value == 'Summary') {

            this.rundata = data;
            this.requestRunDetails(flowcellId);
            requestMinknowMessages(flowcellId);

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
            this.requestCumuBases(flowcellId);
            this.requestHistogramData(flowcellId);
            this.requestQualTime(flowcellId);
            this.requestLengthTime(flowcellId);
            this.requestMaxLengthTime(flowcellId);
            this.requestSpeed(flowcellId);
            this.requestChannelSummaryData(flowcellId);
            this.requestReference(flowcellId);

        } else if (flowcell_selected_tab_input.value == 'Live Event Data') {

            this.rundata = data;
            this.requestRunDetails(flowcellId);
            this.requestLiveRunStats(flowcellId);

        } else {

            this.requestKraken(flowcellId);
            this.requestPafData(flowcellId);
            this.requestGfaData(flowcellId);
        }
    }).bind(this));
};
