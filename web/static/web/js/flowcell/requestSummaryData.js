function requestSummaryData (id) {
    /*
     * Request summary by barcode data
     */
    var url = "/api/v1/flowcells/" + id + "/summarybarcode";

    $.get(url, function (data) {
        if (data.length > 0) {
            var summary = {};
            for (var i = 0; i < data.length; i++) {
                var item = data[i];

                if (summary[item.barcodename] === undefined) {
                    summary[item.barcodename] = {};
                }

                if (summary[item.barcodename][item.typename] === undefined) {
                    summary[item.barcodename][item.typename] = {};
                    summary[item.barcodename][item.typename]["read_count"] = 0;
                    summary[item.barcodename][item.typename]["pass_count"] = 0;
                    summary[item.barcodename][item.typename]["yield"] = 0;
                    summary[item.barcodename][item.typename]["pass_yield"] = 0;
                    summary[item.barcodename][item.typename]["max_length"] = 0;
                    summary[item.barcodename][item.typename]["pass_max_length"] = 0;

                }
                summary[item.barcodename][item.typename]["read_count"] += item.read_count;
                summary[item.barcodename][item.typename]["pass_count"] += item.pass_count;
                summary[item.barcodename][item.typename]["yield"] += item.total_length;
                summary[item.barcodename][item.typename]["pass_yield"] += item.pass_length;
                if (item.max_length > summary[item.barcodename][item.typename]["max_length"]) {
                    summary[item.barcodename][item.typename]["max_length"] = item.max_length;
                }
                if (item.pass_max_length > summary[item.barcodename][item.typename]["pass_max_length"]) {
                    summary[item.barcodename][item.typename]["pass_max_length"] = item.pass_max_length;
                }


            }
            var summaries = {};
            for (var barcodename in summary) {
                if (summaries[barcodename] === undefined) {
                    summaries[barcodename] = {};
                }
                for (var typename in summary[barcodename]) {
                    if (summaries[barcodename][typename] === undefined) {
                        summaries[barcodename][typename] = {
                            "read_count": null,
                            "yield": null,
                            "average_read_length": null,
                            "max_length": null
                        };

                        summaries[barcodename][typename]["read_count"] = {};
                        summaries[barcodename][typename]["read_count"]['all'] = {
                            "name": typename,
                            "data": [summary[barcodename][typename]["read_count"]],//[summaries[item.barcodename][item.typename]["read_count"]["data"] + item.read_count],
                            "animation": false
                        };
                        summaries[barcodename][typename]["read_count"]['pass'] = {
                            "name": typename,
                            "data": [summary[barcodename][typename]["pass_count"]],//[summaries[item.barcodename][item.typename]["read_count"]["data"] + item.read_count],
                            "animation": false
                        };
                        summaries[barcodename][typename]["read_count"]['fail'] = {
                            "name": typename,
                            "data": [(summary[barcodename][typename]["read_count"] - summary[barcodename][typename]["pass_count"])],//[summaries[item.barcodename][item.typename]["read_count"]["data"] + item.read_count],
                            "animation": false
                        };

                        summaries[barcodename][typename]["yield"] = {};
                        summaries[barcodename][typename]["yield"]['all'] = {
                            "name": typename,
                            "data": [summary[barcodename][typename]["yield"]],
                            "animation": false
                        };
                        summaries[barcodename][typename]["yield"]['pass'] = {
                            "name": typename,
                            "data": [summary[barcodename][typename]["pass_yield"]],
                            "animation": false
                        };
                        summaries[barcodename][typename]["yield"]['fail'] = {
                            "name": typename,
                            "data": [(summary[barcodename][typename]["yield"] - summary[barcodename][typename]["pass_yield"])],
                            "animation": false
                        };

                        summaries[barcodename][typename]["average_read_length"] = {};
                        summaries[barcodename][typename]["average_read_length"]['all'] = {
                            "name": typename,
                            "data": [summary[barcodename][typename]["yield"] / summary[barcodename][typename]["read_count"]],
                            "animation": false
                        };
                        summaries[barcodename][typename]["average_read_length"]['pass'] = {
                            "name": typename,
                            "data": [summary[barcodename][typename]["pass_yield"] / summary[barcodename][typename]["pass_count"]],
                            "animation": false
                        };
                        summaries[barcodename][typename]["average_read_length"]['fail'] = {
                            "name": typename,
                            "data": [(summary[barcodename][typename]["yield"] - summary[barcodename][typename]["pass_yield"]) / (summary[barcodename][typename]["read_count"] - summary[barcodename][typename]["pass_count"])],
                            "animation": false
                        };

                        summaries[barcodename][typename]["max_length"] = {};
                        summaries[barcodename][typename]["max_length"]['all'] = {
                            "name": typename,
                            "data": [summary[barcodename][typename]["max_length"]],
                            "animation": false
                        };
                        summaries[barcodename][typename]["max_length"]['pass'] = {
                            "name": typename,
                            "data": [summary[barcodename][typename]["pass_max_length"]],
                            "animation": false
                        };
                        summaries[barcodename][typename]["max_length"]['fail'] = {
                            "name": typename,
                            "data": [summary[barcodename][typename]["max_length"]],
                            "animation": false
                        };

                    }
                }

            }

            self.updateSummaryBasedCharts(summaries);

        }

    });

};