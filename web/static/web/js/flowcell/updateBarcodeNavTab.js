function updateBarcodeNavTab() {
    // Set the request data function to have access to the scope of this function
    var requestData = this.requestData.bind(this);
    // Get the glowcell UD
    var flowcell_id = get_selected_flowcell();
    // get the selected barcode
    var selected_barcode = get_selected_barcode();
    // Get the unorder list element that we have the barcodes under
    var ul = document.getElementById("nav-tabs-barcodes");
    // Remove all the current list elements
    ul.innerHTML = "";
    // get the barcodes
    var sortedBarcodes = this.barcodes;
    // Iterate the array
    for (var i = 0; i < sortedBarcodes.length; i++) {
        // create a list element
        var li = document.createElement("li");
        // Add the barcode-tab class to the new element
        li.classList.add('barcode-tab');
        // create an a element to add to the list element, we're going to attach the click event listener to this tab
        var a = document.createElement("a");
        // add event listener to a element
        a.addEventListener('click', function(event) {
            // Get the text for the selected barcode, upon click
            var selected_barcode = event.target.innerText;
            // set the selected barcode to this newly clicked barcode
            set_selected_barcode(selected_barcode);

            console.log('Calling requestData from inside updateBarcodeNavTab. >>>');
            // call request data, where we will load the data for this barcode for the selected tab
            requestData(flowcell_id);
            console.log('Calling requestData from inside updateBarcodeNavTab. <<<');
            // get all the barcode tabs as elements
            var barcode_tabs = document.querySelectorAll('.barcode-tab');
            // for each of these elements
            barcode_tabs.forEach(function(value) {
                // remove the active class from the tab
                value.classList.remove('active');
                // if it's the selected barcode, add the active class
                if (value.firstChild.text == selected_barcode) {

                    value.classList.add("active");
                }
            });
        });

        // set the a element to the barcode text
        a.text = sortedBarcodes[i];
        // set the barcode class to active, redundant but not terrible
        if (sortedBarcodes[i] == selected_barcode) {
            li.classList.add("active");
        }
        // append the a element to the list element
        li.appendChild(a);
        // append the list element to the unordered list
        ul.appendChild(li);
    }
};
