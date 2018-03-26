function requestData() {

    var flowcellId = this.flowcellId;

    var url_run = '/api/v1/flowcells/' + flowcellId;

    $.get(url_run, (function (data) {

        if (this.summary !== null) {
            if (this.summary["All reads"]['Template']["max_length"]["all"]["data"][0] >= 1000000 && this.millionaire != true) {
                $('#eastermodal').modal('show');
                this.millionaire = true;
            }
        }

        var barcodes = new Set();

        for (var j = 0; j < data.length; j++) {
            for (var i = 0; i < data[j].barcodes.length; i++) {
                barcodes.add(data[j].barcodes[i].name)
            }
        }

        this.barcodes = Array.from(barcodes).sort();
        this.rundata = data;
        this.updateBarcodeNavTab();
        this.requestSummaryData(flowcellId);
        this.requestHistogramData(flowcellId);
        this.requestChannelSummaryData(flowcellId);
        this.requestRunDetails(flowcellId);
        this.requestLiveRunStats(flowcellId);
        this.liveUpdateTasks(flowcellId);
        this.requestKraken(flowcellId);
        this.requestPafData(flowcellId);
        this.requestCumuBases(flowcellId);
        this.requestQualTime(flowcellId);
        this.requestLengthTime(flowcellId);
        this.requestMaxLengthTime(flowcellId);
        this.requestSpeed(flowcellId);
        this.requestGfaData(flowcellId);
        this.requestMessages();

    }).bind(this));

};
