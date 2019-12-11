class MMessageView extends MView {
    // A view of the message - appends a paragraph element to a page - extends the more generic MView class
    // calls the update functionality of the MView to update a HTML element
    constructor(elemento, type = "info") {
        super(elemento);
        this._type = type;
    }

    template(model) {
        // If we have a text element in the model argument return the paragraph element containing the text
        // else return an empty paragraph element
        return model.texto ? `<p class="alert alert-`+ this._type +`">${model.texto}</p>` : '<p></p>';
    }
}
