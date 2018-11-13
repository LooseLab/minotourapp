class FlowcellTabController {

    constructor(flowcell_id) {

        // this._nav_summary_data = document.querySelector('#nav-summary-data');
        // this._nav_live_event_data = document.querySelector('#nav-live-event-data');
        // this._nav_basecalled_data = document.querySelector('#nav-basecalled-data');
        // this._nav_reads_data = document.querySelector('#nav-reads');
        // this._nav_sequence_mapping = document.querySelector('#nav-sequence-mapping');
        // this._nav_sequence_assembly = document.querySelector('#nav-sequence-assembly');
        // this._nav_tab_metagenomics = document.querySelector('#nav-metagenomics');
        // this._nav_tasks = document.querySelector('#nav-tasks');

        this._tabs_view = new MFlowcellTabsView(document.querySelector('.nav-tabs'));

        this._tabs = new FlowcellTabList();

        let tab_summary = new FlowcellTab('Summary', 'summary-data');
        let tab_basecalled_data = new FlowcellTab('Basecalled Data', 'basecalled-data');

        this._tabs.add(tab_summary);
        this._tabs.add(tab_basecalled_data);

        this._tabs_view.update(this._tabs);
        //this._flowcell_services = new FlowcellService();

        //this._tabs = this._flowcell_services.getFlowcellTabs(flowcell_id);


    }

    toggle_tab_content(tab_content_html_element) {

        console.log(`should show ${tab_content_html_element}`);
    }
}



