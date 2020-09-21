class FlowcellTabController {
  /*
    Once a flowcell has been selected, this class is constructed to take charge of the tab switching
     */
  constructor (flowcellId) {
    this._flowcellId = flowcellId
    this._flowcellServices = new FlowcellService()
    this._drawTabs()
    this.tabs = []
    this._first = true
    this._createElementLookupObjects()
    this._redraw_interval = setInterval(() => this._drawTabs(), 60000)
  }

  /**
   * Create a object that has all the tabs and nav HTML objects keyed to the elements ID attribute
   */
  _createElementLookupObjects () {
    const that = this
    this.lookupElement = {}
    $(`.flowcell-tab-ctl`).each((i, obj) => {
      that.lookupElement[obj.id] = obj
    })
  }

  /**
   * Draw the tabs for the first time on the flowcell index.html page, then show one tab of data
   */
  _drawTabs () {
    const promise = this._flowcellServices.getFlowcellTabs(this._flowcellId)
    promise.then((tabs) => {
      const seshTab = getSelectedTab()
      // We are getting the previous tab here and checking it still is available. If it is we show it.
      this.showTabs(tabs)
      this.tabs = tabs
      if (this._first) {
        const activeTab = seshTab !== null && tabs.includes(seshTab) ? seshTab : `summary-data`
        this.toggleTabContent(activeTab)
        this._first = false
      }
    })
  }

  /**
   * Toggle the content to show the tab that the user has selected.
   * @param name {string} The name of the seleceted tab.
   */
  toggleTabContent (name) {
    const navName = `nav-${name}`
    const tabName = `tab-${name}`
    let oldNavTab, oldTabContent, controllerName
    const lookupController = {
      "sequence-mapping": `alignmentController`,
      "basecalled-data": `baseCalledDataController`,
      artic: `articController`,
      tasks: `tasksController`,
      metagenomics: `metagenomicsController`,
      "summary-data": `summaryDataController`,
      "live-event-data": `liveEventController`
    }

    if (getSelectedTab()) {
      oldNavTab = `nav-${getSelectedTab()}`
      oldTabContent = `tab-${getSelectedTab()}`
      // remove active class from current nav tab child element, the link
      // Add hidden tab to old tab content
      this.lookupElement[oldTabContent].classList.add(`hidden`)
      this.lookupElement[oldNavTab].firstElementChild.classList.remove(`active`)
    }
    // set the new tab in Session Storage.
    setSelectedTab(name)
    // add active to newly selected nav time child element, the link
    this.lookupElement[navName].firstElementChild.classList.add(`active`)
    // show the content of the newly selected Div
    this.toggleContent(this.lookupElement[navName], this.lookupElement[tabName])
    controllerName = lookupController[name]
    if (Object.prototype.hasOwnProperty.call(this, controllerName)) {
      // update the tab content by calling the classes update tab method
      console.log(`Updating tab ${name}`)
      this[controllerName].updateTab()
    }
  }

  /**
   * Show the content of the select nav-tab in the tab content div.
   * @param nav {object} HTML object representing the nav tab to be labelled as active.
   * @param tab {object} HTML object representing the tab content to be shown.
   */
  toggleContent (nav, tab) {
    // As of bootstrap 4 we need to add the active class to the a link not the Li element
    nav.children[0].classList.add(`active`)
    tab.classList.remove(`hidden`)
    tab.classList.add(`active`)
    setSelectedTab(tab.id.substr(4))
    // Bad call to request data
    app.requestData(this._flowcellId)
  }

  /**
   * Show the nav tabs on the flowcell_index.html page by removing the hidden class on them.
   * @param tabs {string[]} Array of tab names
   */
  showTabs (tabs) {
    const controllers = {
      "basecalled-data": [`baseCalledDataController`, BasecalledDataController],
      artic: [`articController`, ArticController],
      "sequence-mapping": [`alignmentController`, AlignmentController],
      tasks: [`tasksController`, TasksController],
      metagenomics: [`metagenomicsController`, MetagenomicsController],
      "summary-data": [`summaryDataController`, SummaryTabController],
      "live-event-data": [`liveEventController`, LiveMinKnowController]
    }
    // if we have a controller for this tab, but we are no longer showing it as underlying data has been deleted
    Object.keys(controllers).forEach(controllerName => {
      console.log()
      if (!tabs.includes(controllerName) && this[controllers[controllerName][0]]) {
        // delete the controller
        this.lookupElement[`nav-${controllerName}`].classList.add(`hidden`)
        delete this[controllers[controllerName][0]]
      }
    })
    // Tabs that are shown have their controller initialised here
    tabs.forEach((name) => {
      this.lookupElement[`nav-${name}`].classList.remove(`hidden`)
      // if the name of the tab is in our list of controllers to initialise
      if (Object.keys(controllers).includes(name)) {
        // if we don't have a controller already initialised for this tab
        if (!this[controllers[name][0]]) {
          this[controllers[name][0]] = new controllers[name][1](this._flowcellId)
        }
      }
    })
  }
}
