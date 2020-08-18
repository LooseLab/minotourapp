function addMetaBarcodeTabs (flowcellId, barcodes, tabs) {
  // bind the scope of this function to a the request data function
  var requestData = this.requestData.bind(this)
  // Get the currently selected barcode
  var selected_barcode = getSelectedBarcode(`Metagenomics`)
  // Get the unordered list element
  var ul = document.getElementById(`nav-tabs-meta-barcodes`)
  // Set the unordered list to empty
  ul.innerHTML = ``
  // Get the barcodes, which is an array of barcodes
  var sortedBarcodes = barcodes
  // loop throught the barcodes
  for (var i = 0; i < sortedBarcodes.length; i++) {
    // A lookup dict to decide what colour the barcode tab should be
    var alertLevels = { 0: `green-alert-tab`, 1: `green-alert-tab`, 2: `orange-alert-tab`, 3: `red-alert-tab` }
    // create a list element
    var li = document.createElement(`li`)
    // Get the alert level, 0-3 for each barcode
    var tabLevel = tabs[sortedBarcodes[i]]
    // Get the class for the css styling for the barcode
    var tabLevelClass = alertLevels[tabLevel]
    // Add the class to the tab identifying this as a metagenomics barcode tab
    li.classList.add(`barcode-meta-tab`)
    // If the alert level class is not undefined add it to the list element so it's styled
    if (tabLevelClass !== undefined) {
      li.classList.add(tabLevelClass)
    }
    // create a link element to add to the window
    var a = document.createElement(`a`)
    // add a click event listener to the a element
    a.addEventListener(`click`, function (event) {
      // Get the barcode we have selected from the tabs inner text
      var selected_barcode = event.target.innerText
      // Set the selected barcode as that barcode
      setSelectedBarcode(selected_barcode, `Metagenomics`)
      // Call request data for this flowcell, which will replace all the graphs with data for that barcode
      requestData(flowcellId)
      // Get all the tabs into an array
      var barcode_tabs = document.querySelectorAll(`.barcode-tab`)
      // For each tab loop through, if it's the tab we have selected, set it to active
      barcode_tabs.forEach(function (value) {
        value.classList.remove(`active`)

        if (value.firstChild.text == selected_barcode) {
          value.classList.add(`active`)
        }
      })
    })
    // outside the event listener, so creating original list
    // set the text for the link in the a element to the barcode in this for loop iteration
    a.text = sortedBarcodes[i]
    // if this barcode is the selected barcode, make it active
    if (sortedBarcodes[i] == selected_barcode) {
      li.classList.add(`active`)
    }
    // If this is a barcoded run, otherwise ALl reads is the only barcode in the barcode list we need to display
    if (sortedBarcodes[i] !== `No`) {
      // Add the href link to the list element
      li.appendChild(a)
      // Add the list element to the unordered list
      ul.appendChild(li)
    }
  }
}
