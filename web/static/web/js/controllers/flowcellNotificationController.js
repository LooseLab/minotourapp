class FlowcellNotificationController {
    /**
     * Summary. Controller to handle the display of sent messages and the creation of conditions for messages
     * @class flowcellID
     */
    constructor(flowcellId) {
        this._conditions = {};
        this._flowcellId = flowcellId;
    }


    /**
     * @function Submit our notification choices to the server to create tasks to create objects to pick them up
     */
    submitNotificationChoices(event, data) {

        event.preventDefault();

        // Loop through the input fields and create a dictionary of conditions
        for (const property in Object.keys($("input", data))){
            // Get the input element
            let input = $("input", data)[property];
            console.log(input);
            if (input === undefined){
                continue;
            }
            // If it's a checkbox;
            if (input.type === "checkbox"){
                if (input.checked === true){
                    if (!(input.name in this._conditions)){
                        this._conditions[input.name] = true;
                    }
                }
            }
            // Else if it's the coverage amount we are aiming for
            if (input.name === "coverageAmount"){
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
        }).then(function(response) {
            console.log("sucess");
        }).catch(function(error) {

        });
    }
}