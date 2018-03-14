
function requestData() {

    var flowcellId = this.flowcellId;

    var url_run = '/api/v1/flowcells/' + flowcellId;

    $.get(url_run, (function (data) {

        /*if (self.summary !== null) {
            if (self.summary["All reads"]['Template']["max_length"]["all"]["data"][0] >= 1000000 && self.millionaire != true) {
                $('#eastermodal').modal('show');
                self.millionaire = true;
            }
        }*/

        var barcodes = new Set();

        for (var j = 0; j < data.length; j++) {
            for (var i = 0; i < data[j].barcodes.length; i++) {
                barcodes.add(data[j].barcodes[i].name)
            }
        }

        this.barcodes = Array.from(barcodes).sort();
        this.rundata = data;
        this.updateBarcodeNavTab();
        //self.requestSummaryData(flowcell_id);
        //self.requestHistogramData(flowcell_id);
        //self.requestChannelSummaryData(flowcell_id);
        //self.requestRunDetails(flowcell_id);
        //self.requestLiveRunStats(flowcell_id);
        //self.liveUpdateTasks(flowcell_id);
        //self.requestKraken(flowcell_id);
        //self.requestPafData(flowcell_id);
        //self.requestCumuBases(flowcell_id);
        //self.requestQualTime(flowcell_id);
        //self.requestLengthTime(flowcell_id);
        //self.requestMaxLengthTime(flowcell_id);
        this.requestSpeed(flowcellId);
        //self.requestGfaData(flowcell_id);
        //self.requestMessages();

    }).bind(this));

};

function requestData22(flowcell_id) {

    var url_run = '/api/v1/flowcells/' + flowcell_id;

    $.get(url_run, function (data) {

        /*if (self.summary !== null) {
            if (self.summary["All reads"]['Template']["max_length"]["all"]["data"][0] >= 1000000 && self.millionaire != true) {
                $('#eastermodal').modal('show');
                self.millionaire = true;
            }
        }*/

        var barcodes = new Set();

        for (var j = 0; j < data.length; j++) {
            for (var i = 0; i < data[j].barcodes.length; i++) {
                barcodes.add(data[j].barcodes[i].name)
            }
        }

        this.barcodes = Array.from(barcodes).sort();
        this.rundata = data;
        this.updateBarcodeNavTab().bind(this);
        //self.requestSummaryData(flowcell_id);
        //self.requestHistogramData(flowcell_id);
        //self.requestChannelSummaryData(flowcell_id);
        //self.requestRunDetails(flowcell_id);
        //self.requestLiveRunStats(flowcell_id);
        //self.liveUpdateTasks(flowcell_id);
        //self.requestKraken(flowcell_id);
        //self.requestPafData(flowcell_id);
        //self.requestCumuBases(flowcell_id);
        //self.requestQualTime(flowcell_id);
        //self.requestLengthTime(flowcell_id);
        //self.requestMaxLengthTime(flowcell_id);
        this.requestSpeed(flowcell_id).bind(this);
        //self.requestGfaData(flowcell_id);
        //self.requestMessages();
    });

    console.log(self);

};