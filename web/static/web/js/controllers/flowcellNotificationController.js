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
        this._autoOptionsSelected = {};
        this.options = {
            'backdrop': 'true'
        };
        this._voltageOccuChosen = {};
        this._table = null;
        this.getNotificationsForTable(flowcellId);
    }

    /**
     * @function Get the range choices from the range modal for Occuoancy and Voltage.
     * @param event
     * @param data
     * @private
     */
    _getRangeChoices(event, data) {
        let self = this;
        event.preventDefault();
        const lowerLimit = data[0].value;
        const upperLimit = data[1].value;
        const modal = $("#inputModal");
        console.log(data);
        if (modal.hasClass("Voltage")) {
            self._voltageOccuChosen["Voltage"] = {"lower": lowerLimit, "upper": upperLimit}
        } else if (modal.hasClass("Occupancy")) {
            self._voltageOccuChosen["Occupancy"] = {"lower": lowerLimit, "upper": upperLimit}
        }
        modal.modal("hide")
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
     * @function Disable or enable the fields on the Coverage part of the notification form, based on whether the checkbox is ticked or not.
     * @param event HTML event object
     * @param data The data recovered form the checkbox
     */
    disableCoverageFields(event, data) {
        event.preventDefault();
        // console.log(event);
        // console.log(this);
        if (document.getElementById("coverage").checked) {
            document.getElementById("coverageAmount").disabled = false;
            document.getElementById("reference").disabled = false;
        } else {
            document.getElementById("coverageAmount").disabled = true;
            document.getElementById("reference").disabled = true;
            document.getElementById("contig").disabled = true;
        }
    }

    /**
     * @function Submit our notification choices to the server to create tasks to create objects to pick them up
     */
    submitNotificationChoices(event, data) {
        console.log(data);
        event.preventDefault();

        // Loop through the input fields and create a dictionary of conditions
        const multiInput = document.getElementById("barcodeMulti");
        const multiInputAuto = document.getElementById("autoMulti");
        let refSelect = $("#reference");
        let contigSelect = $("#contig");
        let barcSelect = $("#notification_barcodes");
        let selectedAutoNotifications = multiInputAuto.getValues();
        let justElements = [refSelect, contigSelect, barcSelect];

        let self = this;
        // loop through the select boxes
        if (refSelect[0] !== undefined) {
            let refSelectOptions = refSelect[0].selectedOptions;
            let contigSelectOptions = contigSelect[0].selectedOptions;
            let barcSelectOptions = barcSelect[0].selectedOptions;
            let selects = [[refSelectOptions, "reference"], [contigSelectOptions, "chromosome"], [barcSelectOptions, "barcodes"]];
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
                        console.log(self._barcodesShown);
                        console.log(self._barcodesObj);
                        console.log(selected);
                        // console.log("!!!!!!!!!!!!");
                        // Go through the selected options,
                        let chosen_barcodes = multiInput.getValues();
                        chosen_barcodes.forEach(function (barcode) {
                            self._barcodesShown.forEach(function (elementInSet) {
                                // Get id for barcode to create condition
                                if (typeof elementInSet === "object") {
                                    if (elementInSet[0] === barcode) {
                                        // element in set[1] should be the barcode PK
                                        self._optionsSelected[key].push(elementInSet[1]);
                                    }
                                }
                            });
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
        }
        this._autoOptionsSelected = selectedAutoNotifications;

        let csrftoken = getCookie('csrftoken');

        const axiosInstance = axios.create({
            headers: {"X-CSRFToken": csrftoken}
        });

        axiosInstance.post('/api/v1/messages/conditions', {
            flowcell: this._flowcellId,
            conditions: this._conditions,
            coverage_sets: this._optionsSelected,
            autoChoices: this._autoOptionsSelected,
            rangeValues: this._voltageOccuChosen
        }).then(function (response) {
            // reset the form
            self._optionsSelected["barcodes"] = [];
            document.getElementById("post-form-notification-create").reset();
            // check this is not null, which it is if there is no coverage form
            if (multiInput != null){
                multiInput.deleteAllItems();
            }
            multiInputAuto.deleteAllItems();
            // Set select boxes back to base state
            justElements.forEach(function (selectBox, index) {
                if ([0, 1].includes(index)) {
                    selectBox.append('<option value="" disabled selected hidden>Please Select...</option>');
                }
            });
            // Update the notifications table
            $(".extantNotif").DataTable().ajax.reload(null, false);

            // add event handler for success moda
            $('#basicModal').on('show.bs.modal', function (event) {
                console.log("hide");
                const modal = $(this);
                modal.find('.modal-title').text("Success!");
                let stringy = response.data;

                modal.find('.modal-body h3').html(stringy);
                console.log("hidden")

            });
            $('#basicModal').modal(self.options);

        }).catch(function (error) {
            console.log(error);
            $('#basicModal').on('show.bs.modal', function (event) {
                console.log("hide");
                const modal = $(this);
                modal.find('.modal-title').text("Failed!");
                modal.find('.modal-body h3').text("Failed to create conditions.");
                console.log("hidden")

            });
            $('#basicModal').modal(self.options);
            setTimeout(function () {
                $('#myModal').modal('hide');
            }, 5000);

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
            $("#coverageFormPart").css({
                "visibility": "visible",
                "display": "block"
            });
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
            if (error.response.status === 404) {
                $("#coverageFormPart").remove();
            }
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
            // notificationBarcodes.append(`<option>Select All</option>`);
            // make the barcodes select valid
            document.getElementById("notification_barcodesy").disabled = false;
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
                        notificationBarcodes.append(`<option value="${element[0]}" id="${element[1]}"></option>`);
                    });
                }
            }
        });
    }


    /**
     * @function Get all the Notifications we have already created for this flowcell.
     */
    getNotificationsForTable(flowcellId) {
        let table = $(".extantNotif");
        console.log(table);
        table.DataTable({
            language: {
                infoEmpty: "No records available - Got it?",
            },
            ajax: {
                url: "/api/v1/messages/conditions",
                method: "GET",
                data: {flowcellId}
            },
            columns: [
                {"data": "id"},
                {"data": "notification_type"},
                {"data": "upper_limit"},
                {"data": "lower_limit"},
                {"data": "completed"},
                {"data": "ref_name"},
                {"data": "chrom_line_name"},
                {"data": "coverage_target"},
                {
                    "data": null,
                    "orderable": false,
                    "width": "15%",
                    "render": function (data, type, full) {
                        return `<a class="btn" id="delete_${data.id}" onclick=notificationsController.deleteNotification(${data.id})><i class="fa fa-times"></i> Delete </a>`;
                    }
                }
            ]
        });
    }


    /**
     * @function Delete a notification
     * @param notificationID - The PK of the notification to be deleted
     */
    deleteNotification(notificationID) {
        let csrftoken = getCookie('csrftoken');

        const axiosInstance = axios.create({
            headers: {"X-CSRFToken": csrftoken}
        });
        console.log(notificationID);

        axiosInstance.delete("/api/v1/messages/conditions", {params:{ntpk: notificationID}}).then(response => {
            console.log(response);
            $(".extantNotif").DataTable().ajax.reload(null, false);
        }).catch(error => {

        })
    }
}