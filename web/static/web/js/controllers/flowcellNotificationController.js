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
    }


    /**
     * @function Submit our notification choices to the server to create tasks to create objects to pick them up
     */
    submitNotificationChoices(event, data) {
        console.log(data);
        event.preventDefault();

        // Loop through the input fields and create a dictionary of conditions
        let refSelect = $("#reference")[0].selectedOptions;
        let contigSelect = $("#contig")[0].selectedOptions;
        let barcSelect = $("#notification_barcodes")[0].selectedOptions;
        let selects = [refSelect, contigSelect, barcSelect];
        let self = this;
        selects.forEach(function (selected) {
            if (selected.length > 1) {
                console.log(selected);
                console.log("!!!!!!!!!!!!");
                Object.entries(refSelect, function ([key, value]) {
                    console.log(key, value);
                });
            } else {
                self._optionsSelected[1] = selected[0].innerHTML;
            }
        });

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
            if (input.name === "coverageAmount") {
                this._conditions.coverage = input.value;
            }

        }

        let csrftoken = getCookie('csrftoken');

        const axiosInstance = axios.create({
            headers: {"X-CSRFToken": csrftoken}
        });

        axiosInstance.post('/api/v1/messages/create_conditions', {
            flowcell: this._flowcellId,
            conditions: this._conditions
        }).then(function (response) {
        }).catch(function (error) {
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
                console.log(value);
                if (!(self._referencesSet.has(key))) {
                    $("#reference").append(`<option>${key}</option>`);
                    self._referencesSet.add(key);
                    // loop through contig dictionary
                    self._contigsObj[key] = Object.keys(value);
                    Object.entries(value).forEach(function ([key, value]) {
                        if (!(self._contigSet.has(key))) {
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
            notificationContigs.append(`<option>Select All</option>`);
            var optionSelected = $("option:selected", this);
            console.log("optionSelected");
            console.log(optionSelected);
            self._contigsObj[optionSelected[0].innerHTML].forEach(function (element) {
                self._contigsShown = new Set([...self._contigsShown, ...self._contigsObj[optionSelected[0].innerHTML]]);
                notificationContigs.append(`<option>${element}</option>`);
            });

            let referenceSelected = new Set();
            // for each of our selected contigs
            console.log(Object.entries(optionSelected));
            Object.entries(optionSelected).forEach(([key, value]) => {
                if (value.innerHTML !== undefined) {
                    referenceSelected.add(value.innerHTML);
                    console.log(referenceSelected);
                    console.log(value.innerHTML);
                    console.log(self._contigsObj[value.innerHTML]);
                    self._contigsShown.add(...self._contigsObj[value.innerHTML]);
                    // if it's not already in thered
                    self._contigsObj = new Set([...self._contigsShown].filter(x => !referenceSelected.has(x)));
                }


            });
        });

        $("#contig").on('change', function (e) {
            e.preventDefault();
            let notificationBarcodes = $("#notification_barcodes");
            self._barcodesShown.add("Select All");
            var optionSelected = $("option:selected", this);
            console.log("optionSelected");
            console.log(optionSelected);
            if (optionSelected[0].innerHTML === "Select All") {
                console.log("conditionsal");
                optionSelected = $("#contig")[0].selectedOptions;
                console.log(optionSelected);
            }
            console.log(optionSelected);
            // console.log(Object.entries(optionSelected));

            for (let value of optionSelected) {
                console.log(value);
                if (!value.innerHTML === "Select All") {
                    console.log("No seletc all");
                    continue;
                    console.log(self._barcodesObj);
                    console.log(value);
                    console.log(self._barcodesObj[value.innerHTML]);
                    self._barcodesObj[value.innerHTML].forEach(function (element) {
                        self._barcodesShown = new Set([...self._barcodesShown, ...self._barcodesObj[value.innerHTML]]);
                        notificationBarcodes.append(`<option>${element}</option>`);
                    });
                }

            }


            let contigsSelected = new Set();
            // for each of our selected contigs
            console.log(Object.entries(optionSelected));
            Object.entries(optionSelected).forEach(([key, value]) => {
                if (value.innerHTML !== undefined && value.innerHTML !== "Select All") {
                    contigsSelected.add(value.innerHTML);
                    console.log(contigsSelected);
                    console.log(value.innerHTML);
                    console.log(self._barcodesObj[value.innerHTML]);
                    self._barcodesShown.add(...self._barcodesObj[value.innerHTML]);
                    // if it's not already in thered
                    self._barcodesShown = new Set([...self._barcodesShown].filter(x => !contigsSelected.has(x)));
                }


            });
        });
    }
}