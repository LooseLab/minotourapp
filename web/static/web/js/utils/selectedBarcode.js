/**
 * Return the barcode name for given tab from session storage.
 * @param tab {string} Tab name
 * @returns {string} Barcode name for this tab.
 */
function getSelectedBarcode(tab) {

    return sessionStorage.getItem(`${tab}-barcode`)
}

/**
 * Set the barcode of the tab the user is on in the session storage.
 * @param barcode {string} Barcode name
 * @param tab {string} Tab name
 */
function setSelectedBarcode(barcode, tab) {
    sessionStorage.setItem(`${tab}-barcode`, barcode);
}
