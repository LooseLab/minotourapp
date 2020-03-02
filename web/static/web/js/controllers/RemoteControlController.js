/**
 * Class to handle Ajax requests for the Vue remote control page. Helps to keep things a little more modularised.
 */
class RemoteControlController {

    constructor(vueInstance) {
        console.log("RemoteControlController");
        this.vueInstance = vueInstance;
        this.activeMinions = [];
    }

    /**
     * @function merge the Run stats results into the active minion rsults object held in vueInstance.devices
     * @param minionKey: string of the minion device name.
     * @param runStatsObj: The run stats for this active Minion.
     */
    // _mergeDeviceStats(minionKey, runStatsObj) {
    //     console.log("mErging devices");
    //     let objectIndex = this.vueInstance.devices.findIndex(el => el.minION_name === minionKey);
    //     this.vueInstance.devices[objectIndex] = {...this.vueInstance.devices[objectIndex], ...runStatsObj}
    // }

    /**
     * @function
     * Get active minions for this user. Add them to the data holder of the vue instance on that page.
     */
    getActiveMinionsAndStats() {
        let unique = [];
        let csrftoken = getCookie('csrftoken');
        const axiosInstance = axios.create({
            headers: {"X-CSRFToken": csrftoken}
        });
        let self = this;
        axiosInstance.get("/api/v1/active_minions").then(
            result => {
                console.log(result);
                result.data.forEach(minion => {
                    console.log(minion);
                    self.activeMinions.push(minion.name);
                    this.vueInstance.devices.some(e => e.minION_name === minion.name) ? console.log("hello") : this.vueInstance.devices.push(minion);
                    this.vueInstance.computers.some(e => e.computer === minion.computer) ? console.log("bANtEr") : this.vueInstance.computers.push((({computer, space_available, space_till_shutdown, minknow_version}) => ({computer, space_available, space_till_shutdown, minknow_version}))(minion));
                });
                unique = [...new Set(this.vueInstance.devices.map(item => item.computer))];
                console.log(this.vueInstance.devices);
                // this.getActiveMinionsAndStats()
                // this.getLiveRunInfo(self.activeMinions);
                this.vueInstance.computer_info.number = unique.length
            }
        ).catch(
            error => {
                console.log(error)
            }
        );
    }

    // getLiveRunInfo(activeMinions){
    //     let csrftoken = getCookie('csrftoken');
    //     const axiosInstance = axios.create({
    //         headers: {"X-CSRFToken": csrftoken}
    //     });
    //     console.log(activeMinions);
    //     axiosInstance.get("/api/v1/active_minions/run_stats", {params:{activeMinions: JSON.stringify(activeMinions)}}).then(
    //         result => {
    //             console.log(result);
    //             Object.entries(result.data).forEach(([key, value]) => {
    //                 console.log(key, value);
    //                 this.vueInstance.devices.some(e => e.minION_name === key) ? this._mergeDeviceStats(key, value) : console.log("Hello")
    //             });
    //             console.log(this.vueInstance.devices);
    //
    //         }
    //     ).catch(
    //         error => {
    //             console.log(error)
    //         }
    //     );
// }

}
