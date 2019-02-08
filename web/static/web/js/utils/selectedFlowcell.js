function get_selected_flowcell() {

    var input = document.querySelector('#flowcell-id');
    return input.value;
}

function set_selected_flowcell(flowcell_id) {

    var input = document.querySelector('#flowcell-id');
    input.value = flowcell_id;
}
