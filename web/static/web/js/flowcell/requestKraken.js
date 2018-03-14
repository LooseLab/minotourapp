function requestKraken(id) {
    var parsedkraken = '/api/v1/flowcells/' + id + '/krakenparse/';
    $.get(parsedkraken, function (data) {
        //console.log(data);
        krakendata = [];
        var list = $("#krakenSelectBarcode");
        var list2 = $("#krakenSelectRead");
        for (var i = 0; i < data.length; i++) {
            var barcodename = data[i]["barcode_name"];
            var readtype = data[i]["type_name"]
            if (!(readtype in krakendata)) {
                krakendata[readtype] = [];
                if ($("#krakenSelectRead option[value='" + data[i]["type_name"] + "']").val() === undefined) {
                    list2.append(new Option(data[i]["type_name"], data[i]["type_name"]))
                }
            }
            if (!(barcodename in krakendata[readtype])) {
                krakendata[readtype][barcodename] = [];
                if ($("#krakenSelectBarcode option[value='" + data[i]["barcode_name"] + "']").val() === undefined) {
                    list.append(new Option(data[i]["barcode_name"], data[i]["barcode_name"]));
                }

            }
            //if (data[i]['percentage'] >= 0.05 && data[i]["parent"] != "Input") {
            //if (data[i]['indentation']>8){
            krakendata[readtype][barcodename].push([data[i]["parent"], data[i]["sci_name"], data[i]['percentage'], data[i]['indentation']])
            //}
        }
        //console.log(krakendata);
        list.change(function () {
            var selectedType = list2.find(":selected").text();
            var selectedBarcode = list.find(":selected").text();
            var minimum = $("#lowerboundselect").val();
            //console.log(minimum);
            var slimmedkraken = [];
            for (var i = 0; i < krakendata[selectedType][selectedBarcode].length; i++) {
                if (krakendata[selectedType][selectedBarcode][i][3] >= minimum) {
                    slimmedkraken.push([krakendata[selectedType][selectedBarcode][i][0], krakendata[selectedType][selectedBarcode][i][1], krakendata[selectedType][selectedBarcode][i][2]]);
                }
            }
            ;

            //console.log(selectedBarcode);
            var chart = Highcharts.chart('kraken-sankey', {
                title: {
                    text: 'Kraken ' + selectedType + ' ' + selectedBarcode
                },

                series: [{
                    keys: ['from', 'to', 'weight'],
                    data: slimmedkraken,
                    //data: krakendata[selectedType][selectedBarcode],
                    type: 'sankey',
                    curveFactor: 0,
                    name: 'Kraken Output'
                }]
            });
        });
    })
};
