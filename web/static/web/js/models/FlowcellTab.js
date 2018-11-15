class FlowcellTab {

    constructor(text, name, is_active) {

        this._text = text;
        this._name = name;
        this._is_active = is_active;
    }

    get text() {

        return this._text;
    }

    get name() {

        return this._name;
    }

    get is_active() {

        return this._is_active;
    }

    get html_element() {

        return 'nav-' + this.name;
    }
}