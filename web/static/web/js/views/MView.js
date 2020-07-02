class MView {
  /*
    Create the view with the element as a an argument - the element is a HTML element
     */
  constructor (element) {
    this._element = element
  }

  update (model) {
    // set the element inner html to the passed html#
    const garboCode = $(this._element)
    garboCode.css(`display`, `block`)
    garboCode.html(this.template(model))
    garboCode.fadeIn()
    setTimeout(() => {
      garboCode.fadeOut()
      garboCode.css(`display`, `none`)
    }, 10000)
  }
}
