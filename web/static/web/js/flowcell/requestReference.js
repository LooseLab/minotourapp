function requestReference() {
    var url = "/api/v1/reference/";
    $.get(url, function (data) {
        var references = [];
        for (var i = 0; i < data.length; i++) {
            references.push(data[i]);
        }
        return references;
    });
}
