function requestLiveRunStats(id) {
    //console.log('lastread ' + this.lastread);
    var url_livestats = '/api/v1/flowcells/' + id + '/runstats/' + this.lastread;
    $.get(url_livestats, function (data) {

        if (data.length > 0) {
            self.needtoupdatecharts = true;
            self.lastread = data[data.length - 1].id;
            self.lasttime = new Date(data[data.length - 1].sample_time)
            //console.log(self.rundata.start_time);
            //console.log(self.lasttime.toISOString());

            for (var i = 0; i < data.length; i++) {
                //console.log(data[i]);
                timestamp = new Date(data[i].sample_time).getTime();
                self.livedata.live_read_count = data[i].minKNOW_read_count;
                self.livedata.voltage.push([timestamp, data[i].voltage_value]);
                self.livedata.asictemp.push([timestamp, data[i].asic_temp]);
                self.livedata.heatsinktemp.push([timestamp, data[i].heat_sink_temp]);
                self.livedata.strand.push([timestamp, data[i].strand]);
                self.livedata.good_single.push([timestamp, data[i].good_single]);
                self.livedata.currpercentage = data[i].occupancy;
                self.livedata.currstrand = data[i].strand;
                self.livedata.percentage.push([timestamp, data[i].occupancy]);
                self.livedata.yield_history.push([timestamp, data[i].event_yield]);
                self.livedata.meanratio_history.push([timestamp, data[i].mean_ratio]);
                self.livedata.instrand_history.push([timestamp, data[i].in_strand]);
                self.livedata.openpore_history.push([timestamp, parseInt(data[i].open_pore)]);
                var myStringArray = ["above", "adapter", "below", "good_single", "strand", "inrange", "multiple", "pending_mux_change", "saturated", "unavailable", "unblocking", "unclassified", "unknown"];
                var arrayLength = myStringArray.length;
                //console.log(parseInt(data[i][myStringArray[4]]));
                for (var j = 0; j < arrayLength; j++) {
                    if (isNaN(data[i][myStringArray[j]])) {
                        self.livedata.pore_history[myStringArray[j]].push([timestamp, 0]);
                        //console.log("found a NAN");
                        //console.log(data[i][myStringArray[i]]);
                    } else {
                        self.livedata.pore_history[myStringArray[j]].push([timestamp, parseInt(data[i][myStringArray[j]])]);
                    }

                }
            }
            self.calculatereadtoeventscaling();

            if (self.needtoupdatecharts == true) {
                self.updateLiveHistogram(data);
                self.updateLiveYieldProjection();
                self.updateLiveCumuYield();
                self.updatePoreStats();
                //self.updateTextPredictions();
            }


        }
        //console.log(self.livedata);
    })
};
