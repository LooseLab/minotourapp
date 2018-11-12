class FlowcellTab {

    constructor(text, html_element) {

        this._text = text;
        this._html_element = html_element;
    }

    get text() {

        return this._text;
    }

    get html_element() {

        return this._html_element;
    }
}