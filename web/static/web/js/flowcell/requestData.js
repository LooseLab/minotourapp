
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

        console.log('>> summary');
        console.log(this.summary);

        this.requestSummaryData(flowcellId);

        console.log(this.summary);
        console.log('<< summary');

        this.requestHistogramData(flowcellId);
        this.requestChannelSummaryData(flowcellId);
        this.requestRunDetails(flowcellId);
        this.requestLiveRunStats(flowcellId);
        this.liveUpdateTasks(flowcellId);
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
