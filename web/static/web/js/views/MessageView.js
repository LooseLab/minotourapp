class MessageView extends MView {
  // A view of the message - appends a paragraph element to a page - extends the more generic MView class
  // calls the update functionality of the MView to update a HTML element
  constructor (element, type = `info`) {
    super(element)
    this._type = type
  }

  /**
   *
   * @param model {string}
   * @return {string}
   */
  template (model) {
    // If we have a text element in the model argument return the paragraph element containing the text
    // else return an empty paragraph element
    return model.length ? `<p class="alert alert-` + this._type + `">${model}</p>` : `<p></p>`
  }
}
