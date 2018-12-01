function addMetaBarcodeTabs(flowcellId, barcodes, tabs){
    var requestData = this.requestData.bind(this);

    var selected_barcode = get_selected_barcode();

    var ul = document.getElementById("nav-tabs-meta-barcodes");

    ul.innerHTML = "";

    var sortedBarcodes = barcodes;

    for (var i = 0; i < sortedBarcodes.length; i++) {

        var alertLevels = {0: "green-alert-tab", 1 : "yellow-alert-tab", 2 : "orange-alert-tab", 3 : "red-alert-tab"};

        var li = document.createElement("li");

        var tabLevel = tabs[sortedBarcodes[i]];

        var tabLevelClass = alertLevels[tabLevel];

        li.classList.add('barcode-meta-tab');

        li.classList.add(tabLevelClass);

        var a = document.createElement("a");
        // a.onclick = self.updateChartsBasedOnBarcode;

        a.addEventListener('click', function(event) {

            var selected_barcode = event.target.innerText;

            console.log('clicked on ' + selected_barcode);

            set_selected_barcode(selected_barcode);

            requestData(flowcellId);

            var barcode_tabs = document.querySelectorAll('.barcode-tab');

            barcode_tabs.forEach(function(value) {
                value.classList.remove('active');

                if (value.firstChild.text == selected_barcode) {

                    value.classList.add("active");
                }
            });
        });

        //a.href = "#";
        a.text = sortedBarcodes[i];

        if (sortedBarcodes[i] == selected_barcode) {
            li.classList.add("active");
        }
        li.appendChild(a);
        ul.appendChild(li);
    }
}