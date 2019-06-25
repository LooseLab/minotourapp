class MFlowcellTabsView extends MView {

    constructor(elemento) {
        super(elemento);
    }

    template(model) {

        var response = '';

        response = model._tabs.map(n => {

            if (n.is_active) {

                return '<li role="presentation" class="active"><a class="flowcell-tab" id="' + n.html_element +'" role="tab" onclick="flowcellController.flowcell_tab_controller.toggle_tab_content(\'' + n.name + '\')">' + n.text +'</a></li>'
            } else {

                return '<li role="presentation"><a class="flowcell-tab" id="' + n.html_element +'" role="tab" onclick="flowcellController.flowcell_tab_controller.toggle_tab_content(\'' + n.name + '\')">' + n.text +'</a></li>'
            }
        }).join('')

        return response;
    }
}
