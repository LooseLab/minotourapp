function requestReference() {
    var url = "/api/v1/reference/";
    $.get(url, function (data) {
        var references = [];
        for (var i = 0; i < data.length; i++) {
            //console.log(data[i]);
            references.push(data[i]);
        }
        //console.log(references);
        //this.references = references;
        return references;
    });
}