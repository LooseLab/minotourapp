class FlowcellNotificationController {
    /**
     * Summary. Controller to handle the display of sent messages and the creation of conditions for messages
     * @class flowcellID
     */
    constructor(flowcellId) {
        this._conditions = {};
        this._flowcellId = flowcellId;
        this._referencesSet = new Set();
        this._contigSet = new Set();
        this._contigsObj = {};
        this._contigsShown = new Set();
        this.addHandlerToSelect();
        this._barcodesObj = {};
        this._barcodesShown = new Set();
        this._optionsSelected = {};
        this._addFormClearOnReload();
    }

    /**
     * Clear form when we refresh or reload the page
     * @private
     */
    _addFormClearOnReload() {
        $(window).on("load", function () {
            document.getElementById("post-form-notification-create").reset();
        });
    }


    /**
     * @function Submit our notification choices to the server to create tasks to create objects to pick them up
     */
    submitNotificationChoices(event, data) {
        console.log(data);
        event.preventDefault();

        // Loop through the input fields and create a dictionary of conditions
        let refSelect = $("#reference");
        let contigSelect = $("#contig");
        let barcSelect = $("#notification_barcodes");
        let refSelectOptions = refSelect[0].selectedOptions;
        let contigSelectOptions = contigSelect[0].selectedOptions;
        let barcSelectOptions = barcSelect[0].selectedOptions;
        let justElements = [refSelect, contigSelect, barcSelect];
        let selects = [[refSelectOptions, "reference"], [contigSelectOptions, "chromosome"], [barcSelectOptions, "barcodes"]];
        let self = this;
        // loop through the select boxes

        for (const property in Object.keys($("input", data))) {
            // Get the input element
            let input = $("input", data)[property];
            console.log(input);
            if (input === undefined) {
                continue;
            }
            // If it's a checkbox;
            if (input.type === "checkbox") {
                if (input.checked === true) {
                    if (!(input.name in this._conditions)) {
                        this._conditions[input.name] = null;
                    }
                }
            }
            // Else if it's the coverage amount we are aiming for
            if (input.name === "coverageAmount" && "coverage" in this._conditions) {
                this._conditions.coverage = input.value;
            }

        }
        if ("coverage" in this._conditions) {

            selects.forEach(function (selected) {
                let selections = selected[0];
                // Name of the dict key
                let key = selected[1];
                if (!self._optionsSelected.hasOwnProperty(key)) {
                    self._optionsSelected[key] = [];
                }
                // If it's the barcodes
                if (key === "barcodes") {
                    console.log(selected);
                    console.log("!!!!!!!!!!!!");
                    // Go through the selected options,
                    Object.entries(selections).forEach(function ([_key, value]) {
                        if (value.value !== "Select All"){
                            console.log(key, value.value);
                            self._optionsSelected[key].push(value.value);
                        }

                    });
                } else {
                    console.log(selections);
                    // Check we have the values from the selection drop downs
                    if (selections.length === 0) {
                        throw alert(`Please select a value for ${key}.`);
                    }
                    self._optionsSelected[key] = selections[0].value;
                    console.log(self._optionsSelected);
                }
            });
        }

        let csrftoken = getCookie('csrftoken');

        const axiosInstance = axios.create({
            headers: {"X-CSRFToken": csrftoken}
        });

        axiosInstance.post('/api/v1/messages/create_conditions', {
            flowcell: this._flowcellId,
            conditions: this._conditions,
            coverage_sets: this._optionsSelected
        }).then(function (response) {
            // reset the form
            document.getElementById("post-form-notification-create").reset();
            // Set select boxes back to base state
            justElements.forEach(function (selectBox, index) {
                selectBox.empty();
                if ([0, 1].includes(index)) {
                    selectBox.append('<option value="" disabled selected hidden>Please Select...</option>');
                }
            });
        }).catch(function (error) {
            console.log(error);
            alert("Failed to create conditions.");
        });
    }

    /**
     * @function Get the names and IDs of the references available to the user to display in the autocomplete box
     *
     */
    getReferencesForDataList() {
        let self = this;

        let csrftoken = getCookie('csrftoken');

        const axiosInstance = axios.create({
            headers: {"X-CSRFToken": csrftoken}
        });

        axiosInstance.get(`/api/v1/flowcells/${self._flowcellId}/chromosomecoverage`, {
            params: {names: true}
        }).then(function (names) {
            console.log(names.data);
            // loop through the returned reference names
            Object.entries(names.data).forEach(([key, value]) => {
                console.log(key);
                console.log(value);
                console.log(Object.values(value[0])[0][0]);
                if (!(self._referencesSet.has(key))) {

                    $("#reference").append(`<option value="${value[1]}">${key}</option>`);
                    self._referencesSet.add(key);
                    // ref to contig, with PK as first element to set as value
                    self._contigsObj[key] = [Object.values(value[0])[0][0], ...Object.keys(value[0])];
                    console.log(self._contigsObj);
                    Object.entries(value[0]).forEach(function ([key, value]) {
                        if (!(self._contigSet.has(key))) {
                            console.log(key);
                            console.log(value);
                            self._contigSet.add(key);
                            self._barcodesObj[key] = value;

                        }
                    });
                    // add the barcodes for this contig
                }
            });
        }).catch(function (error) {

        });

    }

    /**
     * @function Add a select all deselect all to the drop down for barcodes on the notifications tab
     */
    addHandlerToSelect() {
        let self = this;
        [$("#notification_barcodes"), $("#contig")].forEach(function (a) {
            a.on("click", function () {
                if ($(this).find(":selected").text() == "Select All") {
                    if ($(this).attr("data-select") == "false")
                        $(this).attr("data-select", "true").find("option").prop("selected", true);
                    else
                        $(this).attr("data-select", "false").find("option").prop("selected", false);
                }
            });
        });

        $("#reference").on('change', function (e) {
            let notificationContigs = $("#contig");
            let bar = $("#notification_barcodes");
            console.log(bar.children());
            console.log(bar.children().length);
            if (bar.children().length > 1) {
                bar.empty();
            }
            notificationContigs.empty();
            notificationContigs.append('<option value="" disabled selected hidden>Please Select...</option>');
            notificationContigs.prop("disabled", false);
            // notificationContigs.append(`<option>Select All</option>`);
            var optionSelected = $("option:selected", this);
            console.log("optionSelected");
            console.log(optionSelected);
            console.log(self._contigsObj);
            console.log(optionSelected[0].innerHTML);
            // self._contigsObj[optionSelected[0].innerHTML].forEach(function (element) {
            //     console.log(element);
            //     self._contigsShown = new Set([...self._contigsShown, ...self._contigsObj[optionSelected[0].innerHTML]]);
            //     // notificationContigs.append(`<option value="${element[0]}">${element[1]}</option>`);
            // });

            let referenceSelected = new Set();
            // for each of our selected contigs
            console.log(Object.entries(optionSelected));
            Object.entries(optionSelected).forEach(([key, value]) => {
                if (value.innerHTML !== undefined) {
                    let a = self._contigsObj[optionSelected[0].innerHTML];
                    console.log(a);
                    self._contigsShown = new Set([...self._contigsShown, ...self._contigsObj[value.innerHTML]]);
                    notificationContigs.append(`<option value="${a[0]}">${a[1]}</option>`);
                    referenceSelected.add(value.innerHTML);
                    console.log(referenceSelected);
                    console.log(value.innerHTML);
                    console.log(self._contigsObj[value.innerHTML]);
                    self._contigsShown.add(...self._contigsObj[value.innerHTML]);
                    // if it's not already in thered
                    self._contigsShown = new Set([...self._contigsShown].filter(x => !referenceSelected.has(x)));
                }


            });
        });

        $("#contig").on('change', function (e) {
            e.preventDefault();
            let notificationBarcodes = $("#notification_barcodes");
            notificationBarcodes.empty();
            // Add select all to the dropdown barcode
            notificationBarcodes.append(`<option>Select All</option>`);
            // Get the selected options of the drop downs
            var optionSelected = $("option:selected", this);
            console.log("optionSelected");
            console.log(optionSelected);
            if (optionSelected[0].innerHTML === "Select All") {
                console.log("conditionsal");
                optionSelected = $("#contig")[0].selectedOptions;
                console.log(optionSelected);
            }
            console.log(optionSelected);

            let contigsSelected = new Set();
            // for each of our selected contigs
            console.log(Object.entries(optionSelected));
            for (let i = 0; i < optionSelected.length; i++) {
                let key = i;
                let value = optionSelected[i];
                console.log(key, value);
                if (value.innerHTML !== undefined && value.innerHTML !== "Select All") {
                    contigsSelected.add(value.innerHTML);
                    console.log(contigsSelected);
                    console.log(value.innerHTML);
                    console.log(self._barcodesObj[value.innerHTML]);
                    self._barcodesShown.add(...self._barcodesObj[value.innerHTML]);
                    // if it's not already in thered
                    self._barcodesShown = new Set([...self._barcodesShown].filter(x => !contigsSelected.has(x)));
                    self._barcodesObj[value.innerHTML].slice(1).forEach(function (element) {
                        self._barcodesShown = new Set([...self._barcodesShown, ...self._barcodesObj[value.innerHTML]]);
                        notificationBarcodes.append(`<option value="${element[1]}">${element[0]}</option>`);
                    });
                }


            }
        });
    }
}