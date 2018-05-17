function write_run_data(textinfo) {

    var text = '';
    text += "<div class='table-responsive table-bordered'>";
    text += "<table class='table'>";
    text += "<tr>";
    text += "<th>Run</th>";
    text += "<th>Run Start Time</th>";
    text += "<th>Last Read Seen</th>";
    text += "<th>MinKNOW Computer Name</th>";
    text += "<th>MinION ID</th>";
    text += "<th>ASIC ID</th>";
    text += "<th>Current Script</th>";
    text += "<th>minKNOW version</th>";
    text += "<th>Flowcell Name</th>";
    text += "<th>Flowcell ID</th>";
    text += "<th>Sample Name</th>";
    text += "<th>Run Name</th>";
    // text += "<th>Sequencing</th>";
    // text += "<th>Barcoded</th>";
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
        text += "<td>" + textinfo[key]["minKNOW_current_script"] + "</td>";
        text += "<td>" + textinfo[key]["minKNOW_version"] + "</td>";
        text += "<td>" + textinfo[key]['name'] + "</td>";
        text += "<td>" + textinfo[key]['flowcellid'] + "</td>";
        text += "<td>" + textinfo[key]["sample_name"] + "</td>";
        text += "<td>" + textinfo[key]['run_name'] + "</td>";
        // if (textinfo[key]['active'] == true) {
        //     text += '<td><i class="fa fa-check" aria-hidden="true"></i></td>';
        // } else {
        //     text += '<td><i class="fa fa-times" aria-hidden="true"></i></td>';
        // }
        // if (textinfo[key]['barcodes'].length > 2) {
        //     text += '<td><i class="fa fa-check" aria-hidden="true"></i></td>';
        // } else {
        //     text += '<td><i class="fa fa-times" aria-hidden="true"></i></td>';
        // }

        text += "</tr>";
    })
    text += "</table></div>";

    $("#minknow-realtime-data-div").html(text);
}


function updatetext(data) {
//<<<<<<< Updated upstream

//    //this.rundetails = data;
    datadump = new Array();

    if (data.length > 0) {
//=======
//    this.rundetails = data;
//    datadump = new Array();
//    console.log(this.rundetails);
//    if (this.rundetails != undefined && data.length > 0) {
//>>>>>>> Stashed changes

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
            // datadump[i].minKNOW_current_script = data[i]['minKNOW_current_script'];
            //datadump[i].name = this.rundata[i]['name'];

            var starttime = new Date(data[i]['minKNOW_start_time']);
            datadump[i].starttime = starttime;
            datadump[i].flowcellid = data[i]['minKNOW_flow_cell_id'];
            //datadump[i].active = this.rundata[i].active;
            //datadump[i].barcodes = this.rundata[i].barcodes;

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
            datadump[i].computer_name = this.rundata[i]['minKNOW_computer'];
            datadump[i].minKNOW_current_script = this.rundata[i]['minKNOW_current_script'];
            var starttime = new Date(Date.parse(this.rundata[i]['start_time']));
            var lastread = new Date(Date.parse(this.rundata[i]['last_read']));
            datadump[i].starttime = starttime;
            datadump[i].lastread = lastread;
            datadump[i].flowcellid = this.rundata[i]['minKNOW_flow_cell_id'];
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
            var div = document.querySelector("#minknow-realtime-data-div");
            div.innerHTML = message;
        }
    }).bind(this))
}
