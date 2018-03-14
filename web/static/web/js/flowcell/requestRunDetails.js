function updatetext() {

    datadump = new Array();

    if (self.rundetails != undefined && self.rundetails.length > 0) {

        for (var i in self.rundetails) {

            if (!self.rundetails.hasOwnProperty(i)) continue;

            datadump[i] = new Array;
            datadump[i].minIONname = self.rundetails[i]['minION_name'];
            datadump[i].minKNOW_asic_id = self.rundetails[i]['minKNOW_asic_id'];
            datadump[i].minKNOW_version = self.rundetails[i]['minKNOW_version'];
            datadump[i].run_id = self.rundetails[i]['minKNOW_hash_run_id'];
            datadump[i].run_name = self.rundetails[i]['minKNOW_run_name'];
            datadump[i].sample_name = self.rundetails[i]['minKNOW_sample_name'];
            datadump[i].computer_name = self.rundetails[i]['minKNOW_computer'];
            datadump[i].minKNOW_current_script = self.rundetails[i]['minKNOW_current_script'];
            datadump[i].name = self.rundata[i]['name'];

            var starttime = new Date(self.rundetails[i]['minKNOW_start_time']);
            datadump[i].starttime = starttime;
            datadump[i].flowcellid = self.rundetails[i]['minKNOW_flow_cell_id'];
            datadump[i].active = self.rundata[i].active;
            datadump[i].barcodes = self.rundata[i].barcodes;

        }

    } else {

        for (var i in self.rundata) {

            datadump[i] = new Array;
            datadump[i].minIONname = self.rundata[i]['minION_name'];
            datadump[i].minKNOW_asic_id = self.rundata[i]['minKNOW_asic_id'];
            datadump[i].minKNOW_version = self.rundata[i]['minKNOW_version'];
            datadump[i].run_id = self.rundata[i]['minKNOW_hash_run_id'];
            datadump[i].run_name = self.rundata[i]['run_name'];
            datadump[i].sample_name = self.rundata[i]['minKNOW_sample_name'];
            datadump[i].computer_name = self.rundata[i]['minKNOW_computer'];
            datadump[i].minKNOW_current_script = self.rundata[i]['minKNOW_current_script'];
            var starttime = new Date(Date.parse(self.rundata[i]['start_time']));
            var lastread = new Date(Date.parse(self.rundata[i]['last_read']));
            datadump[i].starttime = starttime;
            datadump[i].lastread = lastread;
            datadump[i].flowcellid = self.rundata[i]['minKNOW_flow_cell_id'];
            datadump[i].active = self.rundata[i].active;
            datadump[i].barcodes = self.rundata[i].barcodes;

        }

    }

    var sorteddatadump = datadump.sort(dynamicSort("starttime"));
    self.flowcellstart = sorteddatadump[0].starttime;
    self.write_run_data(sorteddatadump);

};

function requestRunDetails(id) {

    var url_RunDetails = '/api/v1/flowcells/' + id + '/rundetails/';

    $.get(url_RunDetails, function (data) {

        self.rundetails = data;

        if (data[0] != undefined) {
            self.livedata.colours_string = data[0].minKNOW_colours_string;
        }

        self.updatetext();

    })

};
