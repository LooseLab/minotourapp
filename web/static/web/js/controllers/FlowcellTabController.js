class FlowcellTabController {
    /*
    Once a flowcell has been selected, this class is constructed to take charge of the tab switchng
     */
    constructor(flowcell_id) {

        this._flowcell_id = flowcell_id;

        this._flowcell_selected_tab = document.querySelector('#flowcell-selected-tab');

        this._nav_summary_data = document.querySelector('#nav-summary-data');
        this._nav_live_event_data = document.querySelector('#nav-live-event-data');
        this._nav_basecalled_data = document.querySelector('#nav-basecalled-data');
        this._nav_reads_data = document.querySelector('#nav-reads');
        this._nav_sequence_mapping = document.querySelector('#nav-sequence-mapping');
        this._nav_advanced_sequence_mapping = document.querySelector('#nav-advanced-sequence-mapping');
        this._nav_sequence_assembly = document.querySelector('#nav-sequence-assembly');
        this._nav_metagenomics = document.querySelector('#nav-metagenomics');
        this._nav_tasks = document.querySelector('#nav-tasks');
        this._nav_sharing = document.querySelector('#nav-sharing');
        this._nav_notifications = document.querySelector('#nav-notifications');
        this._nav_artic = document.querySelector('#nav-artic');

        this._all_navs = [

            this._nav_summary_data,
            this._nav_live_event_data,
            this._nav_basecalled_data,
            this._nav_reads_data,
            this._nav_sequence_mapping,
            this._nav_advanced_sequence_mapping,
            this._nav_sequence_assembly,
            this._nav_metagenomics,
            this._nav_tasks,
            this._nav_sharing,
            this._nav_notifications,
            this._nav_artic
        ];

        this._tab_summary_data = document.querySelector('#tab-summary-data');
        this._tab_live_event_data = document.querySelector('#tab-live-event-data');
        this._tab_basecalled_data = document.querySelector('#tab-basecalled-data');
        this._tab_reads_data = document.querySelector('#tab-reads');
        this._tab_sequence_mapping = document.querySelector('#tab-sequence-mapping');
        this._tab_advanced_sequence_mapping = document.querySelector('#tab-advanced-sequence-mapping');
        this._tab_sequence_assembly = document.querySelector('#tab-sequence-assembly');
        this._tab_metagenomics = document.querySelector('#tab-metagenomics');
        this._tab_tasks = document.querySelector('#tab-tasks');
        this._tab_sharing = document.querySelector('#tab-sharing');
        this._tab_notifications = document.querySelector('#tab-notifications');
        this._tab_artic = document.querySelector('#tab-artic');

        this._all_tabs = [

            this._tab_summary_data,
            this._tab_live_event_data,
            this._tab_basecalled_data,
            this._tab_reads_data,
            this._tab_sequence_mapping,
            this._tab_advanced_sequence_mapping,
            this._tab_sequence_assembly,
            this._tab_metagenomics,
            this._tab_tasks,
            this._tab_sharing,
            this._tab_notifications,
            this._tab_artic
        ];

        this._tabs = new FlowcellTabList();

        this._flowcell_services = new FlowcellService();

        this.draw_tabs();

        this._redraw_interval = setInterval(() => this.redraw_tabs(), 30000);
    }

    draw_tabs(){
        let promise1 = this._flowcell_services.getFlowcellTabs(this._flowcell_id);

        promise1.then((tabs) => {

            this.show_tabs(tabs);
            let seshTab = sessionStorage.getItem("tabName");
            // We are getting the previous tab here and checking it still is available. If it is we show it.
            let activeTab = (seshTab !== null && tabs.includes(seshTab)) ? seshTab : "summary-data";
            this.toggle_tab_content(activeTab);
        });
    }

    redraw_tabs() {

        let promise1 = this._flowcell_services.getFlowcellTabs(this._flowcell_id);

        promise1.then((tabs) => {

            this.show_tabs(tabs);
        });
    }

    toggle_tab_content(name) {

        var nav_name = 'nav-' + name;
        var panel_name = 'panel-' + name;

        this._flowcell_selected_tab.value = nav_name;

        this._all_navs.forEach(nav => {

            nav.children[0].classList.remove('active');
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

            case 'nav-advanced-sequence-mapping':

                this.toggle_content(this._nav_advanced_sequence_mapping, this._tab_advanced_sequence_mapping);
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

            case 'nav-sharing':

                this.toggle_content(this._nav_sharing, this._tab_sharing);
                break;

            case 'nav-notifications':
                this.toggle_content(this._nav_notifications, this._tab_notifications);
                break;

            case 'nav-artic':
                this.toggle_content(this._nav_artic, this._tab_artic);
                break;
        }
    }

    toggle_content(nav, tab) {
        // const topPadding = $(".main-header").height() + 16;
        // $(".content-wrapper").css("padding-top", `${topPadding}px`);
        // As of bootstrap 4 we need to add the active class to the a link not the Li element
        nav.children[0].classList.add('active');
        tab.classList.remove('hidden');
        tab.classList.add('active');
        sessionStorage.setItem("tabName", tab.id.substr(4));
        app.requestData(this._flowcell_id);
    }

    show_tabs(tabs) {

        tabs.forEach((name) => {

            switch(name) {

                case 'summary-data':

                    this._nav_summary_data.classList.remove('hidden');
                    break;

                case 'live-event-data':

                    this._nav_live_event_data.classList.remove('hidden');
                    break;

                case 'basecalled-data':

                    this._nav_basecalled_data.classList.remove('hidden');
                    break;

                case 'reads':

                    this._nav_reads_data.classList.remove('hidden');
                    break;

                case 'sequence-mapping':

                    this._nav_sequence_mapping.classList.remove('hidden');
                    break;

                case 'advanced-sequence-mapping':

                    this._nav_advanced_sequence_mapping.classList.remove('hidden');
                    break;

                case 'sequence-assembly':

                    this._nav_sequence_assembly.classList.remove('hidden');
                    break;

                case 'metagenomics':
                    this._nav_metagenomics.classList.remove('hidden');
                    break;

                case 'tasks':

                    this._nav_tasks.classList.remove('hidden');
                    break;

                case 'sharing':

                    this._nav_sharing.classList.remove('hidden');
                    break;

                case 'notifications':

                    this._nav_notifications.classList.remove("hidden");
                    break;

                case 'artic':

                    this._nav_artic.classList.remove("hidden");
                    break;
            }
        });
    }
}
