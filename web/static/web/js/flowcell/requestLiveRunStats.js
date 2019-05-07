function requestLiveRunStats(id) {

    console.log('This is lastread from flowcellcontroller: ' + flowcell_controller.lastread);

    var url_livestats = '/api/v1/flowcells/' + id + '/runstats/' + flowcell_controller.lastread;

    $.get(url_livestats, (function (result) {

        data = result["data"];
        //console.log(result);
        this.livedata.colours_string = result['minKNOW_colours_string'];
        if (data.length > 0) {

            data.sort(function(a, b){
                var a_date = new Date(a.sample_time);
                var b_date = new Date(b.sample_time);

                return a_date - b_date;
            });

            this.needtoupdatecharts = true;

            flowcell_controller.lastread = data[data.length - 1].id;

            this.lasttime = new Date(data[data.length - 1].sample_time)
            console.log(data);
            for (var i = 0; i < data.length; i++) {
                timestamp = new Date(data[i].sample_time).getTime();
                //this.livedata.colours_string = result['minKNOW_colours_string'];
                this.livedata.live_read_count = data[i].minKNOW_read_count;
                this.livedata.voltage.push([timestamp, data[i].voltage_value]);
                this.livedata.asictemp.push([timestamp, data[i].asic_temp]);
                this.livedata.heatsinktemp.push([timestamp, data[i].heat_sink_temp]);
                this.livedata.strand.push([timestamp, data[i].strand]);
                this.livedata.adapter.push([timestamp,data[i].adapter]);
                this.livedata.pore.push([timestamp,data[i].pore]);
                this.livedata.good_single.push([timestamp, data[i].good_single]);
                this.livedata.currpercentage = data[i].occupancy;
                this.livedata.currstrand = data[i].strand;
                this.livedata.percentage.push([timestamp, data[i].occupancy]);
                this.livedata.yield_history.push([timestamp, data[i].event_yield]);
                this.livedata.meanratio_history.push([timestamp, data[i].mean_ratio]);
                this.livedata.instrand_history.push([timestamp, data[i].in_strand]);
                this.livedata.openpore_history.push([timestamp, parseInt(data[i].open_pore)]);

                var myStringArray = ["above", "adapter", "below", "good_single", "strand", "inrange", "multiple", "pending_mux_change", "saturated", "unavailable", "unblocking", "unclassified", "unknown","zero","pore","no_pore"];
                var arrayLength = myStringArray.length;
                for (var j = 0; j < arrayLength; j++) {
                    //console.log(myStringArray[j]);
                    //console.log(data[i][myStringArray[j]]);
                    if (isNaN(data[i][myStringArray[j]])) {
                        this.livedata.pore_history[myStringArray[j]].push([timestamp, 0]);
                    } else {
                        this.livedata.pore_history[myStringArray[j]].push([timestamp, parseInt(data[i][myStringArray[j]])]);
                    }

                }
            }
            this.calculatereadtoeventscaling();

            if (this.needtoupdatecharts == true) {
                this.updateLiveHistogram(data);
                //this.updateLiveYieldProjection();
                this.updateLiveCumuYield();
                this.updatePoreStats();
                //this.updateTextPredictions();
            }


        }
    }).bind(this));

};
