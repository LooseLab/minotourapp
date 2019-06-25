class MMessage {
    // A message instance - getter and setter for that text
    constructor(texto='') {

        this._texto = texto;
    }

    get texto() {

        return this._texto;
    }

    set texto(texto) {

        this._texto = texto;
    }
}
