class FlowcellTabController {
    /*
    Once a flowcell has been selected, this class is constructed to take charge of the tab switching
     */
    constructor(flowcellId) {
        this._flowcellId = flowcellId;
        this._flowcellServices = new FlowcellService();
        this._drawTabs();
        this.tabs = [];
        this._createElementLookupObjects();
        this._redraw_interval = setInterval(() => this.redrawTabs(), 30000);
    }

    /**
     * Create a object that has all the tabs and nav HTML objects keyed to the elements ID attribute
     */
    _createElementLookupObjects() {
        let that = this
        this.lookupElement = {}
        $(".flowcell-tab-ctl").each((i, obj) => {
            that.lookupElement[obj.id] = obj;
        });
    }

    /**
     * Draw the tabs for the first time on the flowcell index.html page, then show one tab of data
     */
    _drawTabs() {
        let promise = this._flowcellServices.getFlowcellTabs(this._flowcellId);

        promise.then((tabs) => {

            let seshTab = getSelectedTab();
            // We are getting the previous tab here and checking it still is available. If it is we show it.
            let activeTab = seshTab !== null && tabs.includes(seshTab) ? seshTab : "summary-data";
            this.showTabs(tabs);
            this.toggleTabContent(activeTab);
        });
    }

    /**
     * @function Update the tabs on the flowcell index page to include any new tabs sent from the server
     */
    redrawTabs() {
        console.log("redrawing tabs");
        let promise = this._flowcellServices.getFlowcellTabs(this._flowcellId);

        promise.then((tabs) => {
            this.tabs = tabs;
            this.showTabs(tabs);
        });
    }

    /**
     * Toggle the content to show the tab that the user has selected.
     * @param name {string} The name of the seleceted tab.
     */
    toggleTabContent(name) {

        let navName = `nav-${name}`;
        let tabName = `tab-${name}`;
        let oldNavTab, oldTabContent;
        if (getSelectedTab()) {
            oldNavTab = `nav-${getSelectedTab()}`;
            oldTabContent = `tab-${getSelectedTab()}`;
            // remove active class from current nav tab child element, the link
            // Add hidden tab to old tab content
            this.lookupElement[oldTabContent].classList.add("hidden");
            this.lookupElement[oldNavTab].firstElementChild.classList.remove("active");
        }

        // set the new tab in Session Storage.
        setSelectedTab(name);

        // add active to newly selected nav time child element, the link
        this.lookupElement[navName].firstElementChild.classList.add("active");

        // show the content of the newly selected Div
        this.toggleContent(this.lookupElement[navName], this.lookupElement[tabName]);
    }

    /**
     * Show the content of the select nav-tab in the tab content div.
     * @param nav {object} HTML object representing the nav tab to be labelled as active.
     * @param tab {object} HTML object representing the tab content to be shown.
     */
    toggleContent(nav, tab) {
        // As of bootstrap 4 we need to add the active class to the a link not the Li element
        nav.children[0].classList.add("active");
        tab.classList.remove("hidden");
        tab.classList.add("active");
        setSelectedTab(tab.id.substr(4));
        // Bad call to request data
        app.requestData(this._flowcellId);
    }

    /**
     * Show the nav tabs on the flowcell_index.html page by removing the hidden class on them.
     * @param tabs {string[]} Array of tab names
     */
    showTabs(tabs) {
        tabs.forEach((name) => {
            this.lookupElement[`nav-${name}`].classList.remove("hidden");
        });
    }
}
