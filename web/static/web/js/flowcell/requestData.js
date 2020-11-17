function requestData (flowcell_id) {
  console.log(`Requesting data`)
  var flowcellSelectedTabInput = getSelectedTab()
  if (flowcellSelectedTabInput === `reads`) {
    requestReadData(flowcell_id)
  }
};
