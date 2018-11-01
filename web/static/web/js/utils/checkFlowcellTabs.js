function checkFlowcellTabs(flowcell_id) {

    console.log('>>> Inside checkFlowcellTabs');

    var url = "/api/v1/flowcells/" + flowcell_id + "/tabs/";

    var requestData = this.requestData.bind(this);

    // document.querySelector("#nav-"+value.id).addEventListener("click", function(event){
    //     flowcell_selected_tab_input.value = event.target.innerHTML;
    //
    //     requestData(flowcell_id);
    // });


    $.getJSON(url, function (data) {

        console.log(data);

        var items = [];

        $.each(data, function (key, value) {
            items.push(value);
        });

        items.sort(function (a, b) {
            return a.position - b.position;
        });

        var flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');

        // $.each(items, function (key, value) {
        //     /**
        //      if( $(element).length ) returns true when the element is present
        //      if( !$(element).length ) returns true when the element is not present
        //      This is used to prevent duplicate tabs being produced whenever a check is done for new tabs
        //      **/
        //
        //     if (!$('#nav-' + value.id).length) {
        //         var element = '<li><a href="#' + value.id + '"id="nav-' + value.id + '" class="flowcell-tab" data-toggle="tab">' + value.title + '</a></li>';
        //         $('#tasks-li').before(element);
        //         document.querySelector("#nav-"+value.id).addEventListener("click", function(event){
        //             flowcell_selected_tab_input.value = event.target.innerHTML;
        //
        //             requestData(flowcell_id);
        //         });
        //     }
        // });



    });

    setTimeout(checkFlowcellTabs, 60000, flowcell_id);
}

function addStartTabsEvents(flowcellId){

    var requestData = this.requestData.bind(this);
    // add the on click method to the summary and tasks tabs on the page - author Rory
    var flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');

    let tabs = document.querySelectorAll(".flowcell-tab");

    for (let tab of tabs){
        tab.addEventListener("click", function (event) {
            flowcell_selected_tab_input.value = event.target.innerHTML;

            requestData(flowcellId);
        });
    }

}
