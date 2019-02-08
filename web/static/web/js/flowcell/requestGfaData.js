function updateAssemblyCharts(chart, field) {

    while (chart.series.length > 0) {
        chart.series[0].remove();
    }
    for (var barcode of Object.keys(this.assemblySummary)) {
        for (var type of Object.keys(this.assemblySummary[barcode])) {

            chart.addSeries({
                name: barcode + " - " + type,
                data: this.assemblySummary[barcode][type][field]
            });

        }
    }

}

function updateAssemblyBoxplot() {
    var chart = this.ChartBoxPlotContigs;

    var barcats = [];
    var byreadtype = {};

    for (var barcode of Object.keys(this.assemblyLatest)) {

        barcats.push(barcode);

        for (var type of Object.keys(this.assemblyLatest[barcode])) {

            if (byreadtype[type] === undefined) {
                byreadtype[type] = [];
            }
            var contigsizelist = JSON.parse(this.assemblyLatest[barcode][type]['allcontigs']);
            byreadtype[type].push(contigsizelist.sort(function (a, b) {
                return a - b;
            }));
        }
    }

    while (chart.series.length > 0) {
        chart.series[0].remove();
    }

    chart.update({
        xAxis: {
            categories: barcats
        }
    });

    for (var type of Object.keys(byreadtype)) {

        chart.addSeries({
            name: type,
            data: byreadtype[type]
        });

    }
}

function createAssemblyTable() {
    stringtowrite = '<table class="table table-condensed"><tr><th>Barcode</th><th>ReadType</th><th>Input Reads</th><th>Contigs</th><th>Min</th><th>Max</th><th>N50</th><th>Mean</th><th>Total</th><th>Time</th></tr>';
    for (var barcode of Object.keys(this.assemblyLatest)) {
        for (var type of Object.keys(this.assemblyLatest[barcode])) {
            stringtowrite = stringtowrite + '<tr>';
            stringtowrite = stringtowrite + '<td>' + barcode + ' </td>';
            stringtowrite = stringtowrite + '<td>' + type + ' </td>';
            stringtowrite = stringtowrite + '<td>' + this.assemblyLatest[barcode][type]['nreads'] + ' </td>';
            stringtowrite = stringtowrite + '<td>' + this.assemblyLatest[barcode][type]['ncontigs'] + ' </td>';
            stringtowrite = stringtowrite + '<td>' + this.assemblyLatest[barcode][type]['min'] + ' </td>';
            stringtowrite = stringtowrite + '<td>' + this.assemblyLatest[barcode][type]['max'] + ' </td>';
            stringtowrite = stringtowrite + '<td>' + this.assemblyLatest[barcode][type]['n50'] + ' </td>';
            stringtowrite = stringtowrite + '<td>' + this.assemblyLatest[barcode][type]['mean'] + ' </td>';
            stringtowrite = stringtowrite + '<td>' + this.assemblyLatest[barcode][type]['sum'] + ' </td>';
            stringtowrite = stringtowrite + '<td><i>' + new Date(this.assemblyLatest[barcode][type]['time']) + '</i></td> ';
            stringtowrite = stringtowrite + '</tr>';
        }
    }
    stringtowrite = stringtowrite + '</table>';
    document.getElementById('AssemblyTable').innerHTML = stringtowrite;
}

function requestGfaDataCallback(data) {
        if (data.length > 0) {

            var summaries = {}
            var latest = {}

            for (var i = 0; i < data.length; i++) {
                var item = data[i];

                if (summaries[item.barcode_name] === undefined) {
                    summaries[item.barcode_name] = {};
                }

                if (latest[item.barcode_name] === undefined) {
                    latest[item.barcode_name] = {};
                }

                if (summaries[item.barcode_name][item.type_name] === undefined) {
                    summaries[item.barcode_name][item.type_name] = {
                        'ncontigs': [],
                        'n50': [],
                        'sum': []
                    };
                }

                if (latest[item.barcode_name][item.type_name] === undefined) {
                    latest[item.barcode_name][item.type_name] = {};
                }

                summaries[item.barcode_name][item.type_name]['ncontigs'].push([item.nreads, item.ncontigs]);
                summaries[item.barcode_name][item.type_name]['n50'].push([item.nreads, item.n50len]);
                summaries[item.barcode_name][item.type_name]['sum'].push([item.nreads, item.totlen]);

                latest[item.barcode_name][item.type_name]['nreads'] = item.nreads;
                latest[item.barcode_name][item.type_name]['ncontigs'] = item.ncontigs;
                latest[item.barcode_name][item.type_name]['min'] = item.minlen;
                latest[item.barcode_name][item.type_name]['max'] = item.maxlen;
                latest[item.barcode_name][item.type_name]['mean'] = item.meanlen;
                latest[item.barcode_name][item.type_name]['n50'] = item.n50len;
                latest[item.barcode_name][item.type_name]['sum'] = item.totlen;
                latest[item.barcode_name][item.type_name]['time'] = item.timecreated;
                latest[item.barcode_name][item.type_name]['allcontigs'] = item.allcontigs;

            }

            this.assemblySummary = summaries;
            this.assemblyLatest = latest;
            this.updateAssemblyCharts(this.ChartNumContigs, 'ncontigs');
            this.updateAssemblyCharts(this.ChartN50Contigs, 'n50');
            this.updateAssemblyCharts(this.ChartSumContigs, 'sum');

            this.createAssemblyTable();
            this.updateAssemblyBoxplot();

        }

}

function requestGfaData(id) {

    var url = "/api/v1/flowcells/" + id + "/assembly";

    requestGfaDataCallback = requestGfaDataCallback.bind(this);

    $.get(url, requestGfaDataCallback);

}

