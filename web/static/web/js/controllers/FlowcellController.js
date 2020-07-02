class FlowcellController {
/*
This controller is the top level controller, initialised when the flowcell index page is loaded,
it contains references to the other controllers
 */
  constructor (flowcellId) {
    // delete soon
    this._interval = setInterval(requestData, 60000, flowcellId)
    this._flowcellTabController = new FlowcellTabController(flowcellId)
    this._lastRead = 0 // TODO this is a quick fix for the problem on the live data tab. we should refactor this soon.
  }

  get flowcellTabController () {
    return this._flowcellTabController
  }

  get lastRead () {
    return this._lastRead
  }

  set lastRead (val) {
    this._lastRead = val
  }
}
