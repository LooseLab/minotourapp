function updatetext(data) {

    datadump = new Array();

    if (self.rundetails != undefined && data.length > 0) {

        for (var i in data) {

            if (!data.hasOwnProperty(i)) continue;

            datadump[i] = new Array;
            datadump[i].minIONname = data[i]['minION_name'];
            datadump[i].minKNOW_asic_id = data[i]['minKNOW_asic_id'];
            datadump[i].minKNOW_version = data[i]['minKNOW_version'];
            datadump[i].run_id = data[i]['minKNOW_hash_run_id'];
            datadump[i].run_name = data[i]['minKNOW_run_name'];
            datadump[i].sample_name = data[i]['minKNOW_sample_name'];
            datadump[i].computer_name = data[i]['minKNOW_computer'];
            datadump[i].minKNOW_current_script = data[i]['minKNOW_current_script'];
            datadump[i].name = self.rundata[i]['name'];

            var starttime = new Date(data[i]['minKNOW_start_time']);
            datadump[i].starttime = starttime;
            datadump[i].flowcellid = data[i]['minKNOW_flow_cell_id'];
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

    var url = '/api/v1/flowcells/' + id + '/rundetails/';

    $.get(url, function (data) {

        if (data && data.length > 0) {

            self.livedata.colours_string = data[0].minKNOW_colours_string;
            self.updatetext(data);

        }

    })

};
