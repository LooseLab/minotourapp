function write_run_data(textinfo) {

    var text = '';
    text += "<div class='table-responsive'>";
    text += "<table class='table'>";
    text += "<tr>";
    text += "<th>Run</th>";
    text += "<th>Run Start Time</th>";
    text += "<th>Last Read Seen</th>";
    text += "<th>MinKNOW Computer Name</th>";
    text += "<th>MinION ID</th>";
    text += "<th>ASIC ID</th>";
    text += "<th>Sequencing Kit</th>";
    text += "<th>Purpose</th>";
    text += "<th>minKNOW version</th>";
    text += "<th>Flowcell Type</th>";
    text += "<th>Flowcell ID</th>";
    text += "<th>Sample Name</th>";
    text += "<th>Experiment Name</th>";
    text += "</tr>";

    $.each(textinfo, function (key, val) {
        text += "<tr>";
        text += "<td>";
        text += key;
        text += textinfo[key];
        text += "</td>";
        text += "<td>" + textinfo[key]['starttime'] + "</td>";
        text += "<td>" + textinfo[key]['lastread'] + "</td>";
        text += "<td>" + textinfo[key]["computer_name"] + "</td>";
        text += "<td>" + textinfo[key]["minIONname"] + "</td>";
        text += "<td>" + textinfo[key]["minKNOW_asic_id"] + "</td>";
        text += "<td>" + textinfo[key]["sequencing_kit"] + "</td>";
        text += "<td>" + textinfo[key]["minKNOW_exp_script_purpose"].replace('purpose: ','').replace(/"/g,'') + "</td>";
        text += "<td>" + textinfo[key]["minKNOW_version"] + "</td>";
        text += "<td>" + textinfo[key]['flowcelltype'] + "</td>";
        text += "<td>" + textinfo[key]['flowcellid'] + "</td>";
        text += "<td>" + textinfo[key]["sample_name"] + "</td>";
        text += "<td>" + textinfo[key]['experiment_id'] + "</td>";

        text += "</tr>";
    })

    text += "</table></div>";

    $("#run-summaries-div").html(text);

}


function updatetext(data) {

    datadump = new Array();

    if (data.length > 0) {

        for (var i in data) {

            if (!data.hasOwnProperty(i)) continue;
            
            datadump[i] = new Array;
            datadump[i].minIONname = data[i]['minION_name'];
            datadump[i].minKNOW_asic_id = data[i]['minKNOW_asic_id'];
            datadump[i].minKNOW_version = data[i]['minKNOW_version'];
            datadump[i].run_id = data[i]['minKNOW_hash_run_id'];
            datadump[i].run_name = data[i]['minKNOW_run_name'];
            datadump[i].sample_name = data[i]['minKNOW_sample_name'];
            datadump[i].minKNOW_exp_script_purpose = data[i]['minKNOW_exp_script_purpose'];
            datadump[i].sequencing_kit = data[i]['sequencing_kit'];
            datadump[i].computer_name = data[i]['minKNOW_computer'];
            datadump[i].experiment_id = data[i]['experiment_id'];

            var starttime = new Date(data[i]['minKNOW_start_time']);
            datadump[i].starttime = starttime;
            datadump[i].flowcellid = data[i]['minKNOW_flow_cell_id'];
            datadump[i].flowcelltype = data[i]['flowcell_type'];

        }

    } else {

        // TODO the code below shount rely on this.rundata
        for (var i in self.rundata) {

            datadump[i] = new Array;
            datadump[i].minIONname = this.rundata[i]['minION_name'];
            datadump[i].minKNOW_asic_id = this.rundata[i]['minKNOW_asic_id'];
            datadump[i].minKNOW_version = this.rundata[i]['minKNOW_version'];
            datadump[i].run_id = this.rundata[i]['minKNOW_hash_run_id'];
            datadump[i].run_name = this.rundata[i]['run_name'];
            datadump[i].sample_name = this.rundata[i]['minKNOW_sample_name'];
            datadump[i].sequencing_kit = this.rundata[i]['sequencing_kit'];
            datadump[i].minKNOW_exp_script_purpose = this.rundata[i]['minKNOW_exp_script_purpose'];
            datadump[i].computer_name = this.rundata[i]['minKNOW_computer'];
            datadump[i].experiment_id = this.rundata[i]['experiment_id'];
            datadump[i].minKNOW_current_script = this.rundata[i]['minKNOW_current_script'];
            var starttime = new Date(Date.parse(this.rundata[i]['start_time']));
            var lastread = new Date(Date.parse(this.rundata[i]['last_read']));
            datadump[i].starttime = starttime;
            datadump[i].lastread = lastread;
            datadump[i].flowcellid = this.rundata[i]['minKNOW_flow_cell_id'];
            datadump[i].flowcelltype = this.rundata[i]['flowcell_type'];
            datadump[i].active = this.rundata[i].active;
            datadump[i].barcodes = this.rundata[i].barcodes;

        }

    }

    var sorteddatadump = datadump.sort(dynamicSort("starttime"));
    this.flowcellstart = sorteddatadump[0].starttime;
    write_run_data(sorteddatadump);

};


function requestRunDetails(id) {

    var url = '/api/v1/flowcells/' + id + '/rundetails/';

    $.get(url, (function (data) {

        if (data && data.length > 0) {

            this.livedata.colours_string = data[0].minKNOW_colours_string;
            updatetext(data);

        } else {

            var message = "No information available.";
            var div = document.querySelector("#run-summaries-div");
            div.innerHTML = message;

        }

    }).bind(this))

}
