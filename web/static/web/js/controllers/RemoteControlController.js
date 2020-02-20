/**
 * Class to handle Ajax requests for the Vue remote control page. Helps to keep things a little more modularised.
 */
class RemoteControlController {

    constructor(vueInstance) {
        console.log("RemoteControlController");
        this.vueInstance = vueInstance;
    }


    /**
     * @function
     * Get active minions for this user. Add them to the data holder of the vue instance on that page.
     */
    getActiveMinions() {
        let csrftoken = getCookie('csrftoken');
        const axiosInstance = axios.create({
            headers: {"X-CSRFToken": csrftoken}
        });
        axiosInstance.get("/api/v1/activeminions").then(
            result => {
                console.log(result);
                result.data.forEach(minion => {
                    this.vueInstance.devices.includes(minion.name) ? console.log("hello") : this.vueInstance.devices.push(minion);
                });
                console.log(this.vueInstance.devices);

            }
        ).catch(
            error => {
                console.log(error)
            }
        );
    }

    getLiveRunInfo(){
        let csrftoken = getCookie('csrftoken');
        const axiosInstance = axios.create({
            headers: {"X-CSRFToken": csrftoken}
        });
        axiosInstance.get("/api/v1/").then(
            result => {
                console.log(result);
                result.data.forEach(minion => {
                    this.vueInstance.devices.includes(minion.name) ? console.log("hello") : this.vueInstance.devices.push(minion);
                });
                console.log(this.vueInstance.devices);

            }
        ).catch(
            error => {
                console.log(error)
            }
        );
    }

}
