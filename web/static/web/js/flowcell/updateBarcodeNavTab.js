function updateBarcodeNavTab() {

    var selected_barcode = get_selected_barcode();

    var ul = document.getElementById("nav-tabs-barcodes");

    ul.innerHTML = "";

    var sortedBarcodes = this.barcodes;

    for (var i = 0; i < sortedBarcodes.length; i++) {
        var li = document.createElement("li");
        li.classList.add('barcode-tab');
        var a = document.createElement("a");
        // a.onclick = self.updateChartsBasedOnBarcode;

        a.addEventListener('click', function(event) {

            var selected_barcode = event.target.innerText;
            console.log('clicked on ' + selected_barcode);
            set_selected_barcode(selected_barcode);

            var flowcell_id = get_selected_flowcell();
            requestData(flowcell_id);

            var barcode_tabs = document.querySelectorAll('.barcode-tab');
            barcode_tabs.forEach(function(value) {
                value.classList.remove('active');

                if (value.firstChild.text == selected_barcode) {

                    value.classList.add("active");
                }
            });
        });

        a.href = "#";
        a.text = sortedBarcodes[i];

        if (sortedBarcodes[i] == selected_barcode) {
            li.classList.add("active");
        }

        li.appendChild(a);
        ul.appendChild(li);
    }
};