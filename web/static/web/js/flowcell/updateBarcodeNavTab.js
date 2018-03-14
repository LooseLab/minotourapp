function updateBarcodeNavTab() {
    var ul = document.getElementById("nav-tabs-barcodes");

    ul.innerHTML = "";

    var sortedBarcodes = this.barcodes;

    for (var i = 0; i < sortedBarcodes.length; i++) {
        var li = document.createElement("li");
        var a = document.createElement("a");
        a.onclick = self.updateChartsBasedOnBarcode;
        a.href = "#";
        a.text = sortedBarcodes[i];

        if (sortedBarcodes[i] === self.selectedBarcode) {
            li.classList.add("active");
        }

        li.appendChild(a);
        ul.appendChild(li);
    }
};