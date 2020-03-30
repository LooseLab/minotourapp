// /**
//  * Class to handle Ajax requests for the Vue remote control page. Helps to keep things a little more modularised.
//  */
// class RemoteControlController {
//
//     constructor(vueInstance) {
//         console.log("RemoteControlController");
//         this.vueInstance = vueInstance;
//         this.activeMinions = [];
//         // this._interval = setInterval(this.getActiveMinionsAndStats.bind(this), 1000);
//         // this._clearTimer();
//
//     }
//
//     /**
//      * @function Clear a set timeout object when we navigate away from the page, so we aren't constantly calling getActiveMInionAndStats
//      */
//     _clearTimer() {
//         console.log("clearing interval");
//         $(window).on("unload", () => clearInterval(this._interval));
//     }
//
//
//     /**
//      * @function merge the Run stats results into the active minion results object held in vueInstance.devices
//      * @param minion: Dictionary containing all the information relating to the minion
//      * @param vueInstance: access to the vueInstance instance
//      */
//     _mergeDeviceStats(vueInstance, minion) {
//         let objectIndex = vueInstance.devices.findIndex(el => el.minION_name === minion.minIONname);
//         console.log(vueInstance.devices[objectIndex]);
//         vueInstance.devices[objectIndex] = {...vueInstance.devices[objectIndex], ...minion};
//         console.log(vueInstance.devices[objectIndex])
//     }
//
//     _checkStillActive
//
//     /**
//      * @function
//      * Get active minions for this user. Add them to the data holder of the vue instance on that page.
//      */
//     getActiveMinionsAndStats() {
//         let unique = [];
//         let csrftoken = getCookie('csrftoken');
//         const axiosInstance = axios.create({
//             headers: {"X-CSRFToken": csrftoken}
//         });
//         let self = this;
//         axiosInstance.get("/api/v1/active_minions").then(
//             result => {
//                 result.data.forEach(minion => {
//                     self.activeMinions.push(minion.name);
//                     self.vueInstance.devices.some(e => e.minION_name === minion.name) ? self._mergeDeviceStats(self.vueInstance, minion) : self.vueInstance.devices.push(minion);
//                     self.vueInstance.computers.some(e => e.computer === minion.computer) ? null : self.vueInstance.computers.push(
//                         (
//                             (
//                                 {computer, space_available, space_till_shutdown, minknow_version}) => (
//                                 {
//                                     computer,
//                                     space_available,
//                                     space_till_shutdown,
//                                     minknow_version
//                                 }
//                                 )
//                         )
//                         (minion)
//                     );
//                 });
//                 unique = [...new Set(self.vueInstance.devices.map(item => item.computer))];
//                 self.vueInstance.computer_info.number = unique.length;
//                 // t = setTimeout(this.getActiveMinionsAndStats(), 300000 );
//                 // this._clearTimer(t)
//
//             }
//         ).catch(
//             error => {
//                 console.log(error)
//             }
//         );
//
//
//         // getLiveRunInfo(activeMinions){
//         //     let csrftoken = getCookie('csrftoken');
//         //     const axiosInstance = axios.create({
//         //         headers: {"X-CSRFToken": csrftoken}
//         //     });
//         //     console.log(activeMinions);
//         //     axiosInstance.get("/api/v1/active_minions/run_stats", {params:{activeMinions: JSON.stringify(activeMinions)}}).then(
//         //         result => {
//         //             console.log(result);
//         //             Object.entries(result.data).forEach(([key, value]) => {
//         //                 console.log(key, value);
//         //                 this.vueInstance.devices.some(e => e.minION_name === key) ? this._mergeDeviceStats(key, value) : console.log("Hello")
//         //             });
//         //             console.log(this.vueInstance.devices);
//         //
//         //         }
//         //     ).catch(
//         //         error => {
//         //             console.log(error)
//         //         }
//         //     );
// // }
//
//     }
// }