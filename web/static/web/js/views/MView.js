class MView {
    /*
    Create the view with the element as a an argument - the element is a HTML element
     */
    constructor(elemento) {

        this._elemento = elemento;
    }

    template() {
        // not implemented yet
        throw new Error('O método template deve ser implementado');
    }

    update(model) {
        // set the element inner html to the passed html
        this._elemento.innerHTML = this.template(model);
    }
}
