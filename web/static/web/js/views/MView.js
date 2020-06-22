class MView {
  /*
    Create the view with the element as a an argument - the element is a HTML element
     */
  constructor (elemento) {
    this._elemento = elemento
  }

  update (model) {
    // set the element inner html to the passed html#
    const garboCode = $(this._elemento)
    this._elemento.innerHTML = this.template(model)
    garboCode.fadeIn()
    setTimeout(() => {
      garboCode.fadeOut()
    }, 3000)
  }
}
