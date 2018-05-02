function get_selected_barcode() {

    var input = document.querySelector('#flowcell-selected-barcode');
    return input.value;
}

function set_selected_barcode(barcode) {

    var input = document.querySelector('#flowcell-selected-barcode');
    input.value = barcode;
}
