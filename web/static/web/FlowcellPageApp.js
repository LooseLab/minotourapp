
var FlowcellPageApp = {

  init: function () {
    console.log(`initialising flowcellpageapp`)


    this.requestData = requestData

    var flowcell_id = getSelectedFlowcell()


    // ('>> calling request data');
    console.log(`Calling request data from monitor_app. >>>`)
    this.requestData(flowcell_id)
    console.log(`Calling request data from monitor_app. <<<`)
  } // end of init
}
