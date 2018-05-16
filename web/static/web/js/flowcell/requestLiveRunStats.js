function requestLiveRunStats(id) {
    //console.log('lastread ' + this.lastread);
    var url_livestats = '/api/v1/flowcells/' + id + '/runstats/' + this.lastread;

    $.get(url_livestats, (function (data) {

        if (data.length > 0) {

            data.sort(function(a, b){
                var a_date = new Date(a.sample_time);
                var b_date = new Date(b.sample_time);

                return a_date - b_date;
            });

            this.needtoupdatecharts = true;
            this.lastread = data[data.length - 1].id;
            this.lasttime = new Date(data[data.length - 1].sample_time)
            //console.log(this.rundata.start_time);
            //console.log(this.lasttime.toISOString());

            for (var i = 0; i < data.length; i++) {
                //console.log(data[i]);
                timestamp = new Date(data[i].sample_time).getTime();
                this.livedata.live_read_count = data[i].minKNOW_read_count;
                this.livedata.voltage.push([timestamp, data[i].voltage_value]);
                this.livedata.asictemp.push([timestamp, data[i].asic_temp]);
                this.livedata.heatsinktemp.push([timestamp, data[i].heat_sink_temp]);
                this.livedata.strand.push([timestamp, data[i].strand]);
                this.livedata.good_single.push([timestamp, data[i].good_single]);
                this.livedata.currpercentage = data[i].occupancy;
                this.livedata.currstrand = data[i].strand;
                this.livedata.percentage.push([timestamp, data[i].occupancy]);
                this.livedata.yield_history.push([timestamp, data[i].event_yield]);
                this.livedata.meanratio_history.push([timestamp, data[i].mean_ratio]);
                this.livedata.instrand_history.push([timestamp, data[i].in_strand]);
                this.livedata.openpore_history.push([timestamp, parseInt(data[i].open_pore)]);
                var myStringArray = ["above", "adapter", "below", "good_single", "strand", "inrange", "multiple", "pending_mux_change", "saturated", "unavailable", "unblocking", "unclassified", "unknown"];
                var arrayLength = myStringArray.length;
                //console.log(parseInt(data[i][myStringArray[4]]));
                for (var j = 0; j < arrayLength; j++) {
                    if (isNaN(data[i][myStringArray[j]])) {
                        this.livedata.pore_history[myStringArray[j]].push([timestamp, 0]);
                        //console.log("found a NAN");
                        //console.log(data[i][myStringArray[i]]);
                    } else {
                        this.livedata.pore_history[myStringArray[j]].push([timestamp, parseInt(data[i][myStringArray[j]])]);
                    }

                }
            }
            this.calculatereadtoeventscaling();

            if (this.needtoupdatecharts == true) {
                this.updateLiveHistogram(data);
                this.updateLiveYieldProjection();
                this.updateLiveCumuYield();
                this.updatePoreStats();
                //this.updateTextPredictions();
            }


        }
        //console.log(this.livedata);
    }).bind(this));

};
