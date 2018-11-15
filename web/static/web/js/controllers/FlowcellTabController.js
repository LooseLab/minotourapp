class FlowcellTabController {

    constructor(flowcell_id) {

        this._flowcell_id = flowcell_id;

        this._flowcell_selected_tab = document.querySelector('#flowcell-selected-tab');

        this._nav_summary_data = document.querySelector('#nav-summary-data');
        this._nav_live_event_data = document.querySelector('#nav-live-event-data');
        this._nav_basecalled_data = document.querySelector('#nav-basecalled-data');
        this._nav_reads_data = document.querySelector('#nav-reads');
        this._nav_sequence_mapping = document.querySelector('#nav-sequence-mapping');
        this._nav_sequence_assembly = document.querySelector('#nav-sequence-assembly');
        this._nav_metagenomics = document.querySelector('#nav-metagenomics');
        this._nav_tasks = document.querySelector('#nav-tasks');

        this._all_navs = [

            this._nav_summary_data,
            this._nav_live_event_data,
            this._nav_basecalled_data,
            this._nav_reads_data,
            this._nav_sequence_mapping,
            this._nav_sequence_assembly,
            this._nav_metagenomics,
            this._nav_tasks
        ];

        this._tab_summary_data = document.querySelector('#tab-summary-data');
        this._tab_live_event_data = document.querySelector('#tab-live-event-data');
        this._tab_basecalled_data = document.querySelector('#tab-basecalled-data');
        this._tab_reads_data = document.querySelector('#tab-reads');
        this._tab_sequence_mapping = document.querySelector('#tab-sequence-mapping');
        this._tab_sequence_assembly = document.querySelector('#tab-sequence-assembly');
        this._tab_metagenomics = document.querySelector('#tab-metagenomics');
        this._tab_tasks = document.querySelector('#tab-tasks');

        this._all_tabs = [

            this._tab_summary_data,
            this._tab_live_event_data,
            this._tab_basecalled_data,
            this._tab_reads_data,
            this._tab_sequence_mapping,
            this._tab_sequence_assembly,
            this._tab_metagenomics,
            this._tab_tasks
        ];

        // this._tabs_view = new MFlowcellTabsView(document.querySelector('.flowcell-tabs'));

        this._tabs = new FlowcellTabList();

        // this._tabs_view.update(this._tabs);
        // this._flowcell_services = new FlowcellService();

        // this._tabs = this._flowcell_services.getFlowcellTabs(flowcell_id);

        this.toggle_tab_content('summary-data');
    }

    toggle_tab_content(name) {

        var nav_name = 'nav-' + name;
        var panel_name = 'panel-' + name;

        this._flowcell_selected_tab.value = nav_name;

        this._all_navs.forEach(nav => {

            nav.classList.remove('active');
        });

        this._all_tabs.forEach(tab => {

            tab.classList.add('hidden');
        });

        switch (nav_name) {

            case 'nav-summary-data':

                this.toggle_content(this._nav_summary_data, this._tab_summary_data);
                break;

            case 'nav-live-event-data':

                this.toggle_content(this._nav_live_event_data, this._tab_live_event_data);
                break;

            case 'nav-basecalled-data':

                this.toggle_content(this._nav_basecalled_data, this._tab_basecalled_data);
                break;

            case 'nav-reads':

                this.toggle_content(this._nav_reads_data, this._tab_reads_data);
                break;

            case 'nav-sequence-mapping':

                this.toggle_content(this._nav_sequence_mapping, this._tab_sequence_mapping);
                break;

            case 'nav-sequence-assembly':

                this.toggle_content(this._nav_sequence_assembly, this._tab_sequence_assembly);
                break;

            case 'nav-metagenomics':

                this.toggle_content(this._nav_metagenomics, this._tab_metagenomics);
                break;

            case 'nav-tasks':

                this.toggle_content(this._nav_tasks, this._tab_tasks);
                break;

        }
    }

    toggle_content(nav, tab) {

        nav.classList.add('active');
        tab.classList.remove('hidden');
        tab.classList.add('show');
        app.requestData(this._flowcell_id);
    }
}
