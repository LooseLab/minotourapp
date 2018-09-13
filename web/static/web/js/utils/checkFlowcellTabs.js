function checkFlowcellTabs(flowcell_id) {

    console.log("Checking flowcell tabs.");

    var url = "/api/v1/flowcells/tabs/" + flowcell_id;

    var requestData = this.requestData.bind(this);

    $.getJSON(url, function (data) {

        var items = [];

        $.each(data, function (key, value) {
            items.push(value);
        });

        items.sort(function (a, b) {
            return a.position - b.position;
        });

        $.each(items, function (key, value) {
            /**
             if( $(element).length ) returns true when the element is present
             if( !$(element).length ) returns true when the element is not present
             This is used to prevent duplicate tabs being produced whenever a check is done for new tabs
             **/

            if (!$('#nav-' + value.id).length) {
                var element = '<li><a href="#' + value.id + '"id="nav-' + value.id + '" class="flowcell-tab" data-toggle="tab">' + value.title + '</a></li>';
                $('#tasks-li').before(element);
                element
            }
        });

        var flowcell_selected_tab_input = document.querySelector('#flowcell-selected-tab');

        var tabs = document.querySelectorAll('.flowcell-tab');

        for (var i = 0; i < tabs.length; i++) {
            console.log(tabs[i]);
            tabs[i].addEventListener('click', function (event) {

                flowcell_selected_tab_input.value = event.target.innerHTML;
                console.log('Clicked ' + event.target.innerHTML);
                requestData(flowcell_id);
            });
        }
    });

    // setTimeout(checkFlowcellTabs, 60000, flowcell_id);
}
