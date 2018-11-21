class FlowcellTabList {

    constructor() {

        this._tabs = [];
    }

    add(tab) {

        this._tabs.push(tab);
    }

    get tabs() {

        return [].concat(this._tabs);
    }
}