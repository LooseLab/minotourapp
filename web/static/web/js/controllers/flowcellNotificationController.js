class FlowcellNotificationController {
    /**
     * Summary. Controller to handle the display of sent messages and the creation of conditions for messages
     * @class flowcellID
     */
    constructor(flowcellId) {
        this._conditions = {};
        this._flowcellId = flowcellId;
        this._referencesList = [];
        this.addHandlerToSelect();
        this._barcodesObj = {};
        this._barcodesShown = new Set();
    }


    /**
     * @function Submit our notification choices to the server to create tasks to create objects to pick them up
     */
    submitNotificationChoices(event, data) {

        event.preventDefault();

        // Loop through the input fields and create a dictionary of conditions
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
            // loop through the returned reference names
            Object.entries(names.data).forEach(([key, value]) => {
                if (!(self._referencesList.includes(key))) {
                    $("#reference").append(`<option>${key}</option>`);
                    self._referencesList.push(key);

                    // add the barcodes for this contig
                    self._barcodesObj[key] = value;
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
        $("#notification_barcodes").on("click", function () {
            if ($(this).find(":selected").text() == "Select All") {
                if ($(this).attr("data-select") == "false")
                    $(this).attr("data-select", "true").find("option").prop("selected", true);
                else
                    $(this).attr("data-select", "false").find("option").prop("selected", false);
            }
        });
        $("#reference").on('change', function (e) {
            var optionSelected = $("option:selected", this);
            console.log(optionSelected);
            let contigsSelected = new Set();
            // for each of our selected contigs
            Object.entries(optionSelected).forEach(([key, value]) => {
                contigsSelected.add(value.innerHTML);
                self._barcodesShown.add(value.innerHTML);
                self._barcodesShown = new Set([...self._barcodesShown].filter(x => !contigsSelected.has(x)));


            });
        });
    }
}