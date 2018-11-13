class MFlowcellTabsView extends MView {

    constructor(elemento) {
        super(elemento);
    }

    template(model) {

        return `
            ${model._tabs.map(n => `
                <li>
                    <a class="flowcell-tab" id="nav-${n.html_element}" href="#panel-${n.html_element}" data-toggle="tab" role="tab" aria-controls="panel-${n.html_element}" onclick="flowcell_controller.flowcell_tab_controller.toggle_tab_content('${n.html_element}')">${n.text}</a>
                </li>
            
            `).join('')}
        `;
    }
}