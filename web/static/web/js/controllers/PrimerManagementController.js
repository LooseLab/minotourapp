class PrimerManagementController {
  /**
   * Controller for the mapping tab.
   */
  constructor () {
    this._axiosInstance = axios.create({
      headers: { 'X-CSRFToken': getCookie(`csrftoken`) }
    })
    console.log(`this`)
    this._formData = new FormData()
    // no bigger than 100 Mb files
    this._maxFileSize = 100000000
    this._datatableObj = $(`#scheme-table`)
    this._addEventsToDropZone()
    this._addEventsToSubmit()
    this._getUsersSchemesTable()
  }

  _addEventsToSubmit () {
    $(`#submit-scheme`).on(`click`, (event) => {
      this._submitScheme()
    })
  }

  _deleteScheme () {
    console.log(`todo`)
  }

  _getUsersSchemesTable () {
    if ($.fn.DataTable.isDataTable(this._datatableObj)) {
      this._datatableObj.DataTable().ajax.reload(null, false)
    } else {
      this._datatableObj.DataTable({
        createdRow: (row, data, dataIndex, cells) => {
          $(row).find(`.ref-deletion`).on(`click`, (event) => {
            this._deleteScheme(event.currentTarget.getAttribute(`data-delete`))
          })
        },
        ajax: {
          url: `/api/v1/artic/schemes`,
          async: true,
          error: (xhr, error, code) => {
            console.error(xhr)
            console.error(code)
          }
        },
        columnDefs: [
          {
            targets: 0,
            data: `scheme_species`

          },
          {
            targets: 1,
            data: `scheme_versions`
          },
          {
            targets: 2,
            data: `private`,
            render: (data, type, full, meta) => data.private ? `Private` : `Public`
          },
          {
            width: `1.5rem`,
            sortable: false,
            targets: 3,
            data: null,
            render: (data, type, full, meta) => data.deletable ? `<a class="btn icon-task ref-deletion" data-delete="${data.id}"><i class="fa fa-minus-square"></i> Delete </a>` : ``
          }
        ]
      })
    }
  }

  _dropHandler (ev, dropped) {
    console.log(`hello`)
    ev.preventDefault()
    ev.stopImmediatePropagation()
    // check the md5 checksum here
    const formData = this._formData
    const fileApiSystem = dropped ? ev.dataTransfer.items : ev.currentTarget.files
    $(`#drop-zone`).removeClass(`drop-zone-entered`)
    console.log(fileApiSystem)
    if (fileApiSystem.length) {
      if (fileApiSystem.length > 2) {
        alert(`Please only upload two files, scheme.bed and reference.fasta.`)
        return
      }
      const schemeNames = new Set()
      // Use DataTransferItemList interface to access the file(s)
      for (var i = 0; i < fileApiSystem.length; i++) {
        // If dropped items aren't files, reject them
        if (fileApiSystem[i].kind === `file` || !dropped) {
          const file = dropped ? fileApiSystem[i].getAsFile() : fileApiSystem[i]
          console.log(file)
          if (file.size > this._maxFileSize) {
            alert(`${file.name} is too large (${file.size})! Upload refused.`)
          }
          if (file.name.endsWith(`reference.fasta`)) { $(`#reference-file-name`).html(file.name) }
          if (file.name.endsWith(`scheme.bed`)) { $(`#bed-file-name`).html(file.name) }
          schemeNames.add(file.name.split(`.`)[0])
          formData.append(`file_location`, file, file.name.replace(` `, `_`))
          formData.append(`file_name`, file.name.replace(` `, `_`))
        }
      }
      // check files are called the same thing
      console.log(schemeNames)
      if (schemeNames.size === 1) { $(`#scheme-name`).val(Array.from(schemeNames)[0]) } else { alert(`Files must share same name`) }
      $(`#submit-scheme`).prop(`disabled`, false)
    }
  }

  _validateSchemeFormData (formData) {
    // check that all fields have same name, have data
  }

  _submitScheme () {
    // get the form values for the scheme
    $(`#scheme-settings`).find(`.form-group`).find(`input`).each((index, field) => {
      const value = field.type === `checkbox` ? $(field).prop(`checked`) : $(field).val()
      this._formData.append(field.id.replaceAll(`-`, `_`), value)
    })
    this._axiosInstance({
      url: `/api/v1/artic/schemes/`,
      data: this._formData,
      onUploadProgress: progressEvent => {
        const totalLength = progressEvent.lengthComputable ? progressEvent.total : progressEvent.target.getResponseHeader(`content-length`) || progressEvent.target.getResponseHeader(`x-decompressed-content-length`)
        if (totalLength !== null) {
          $(`#uploadProgressBar`).css(`width`, `${Math.round((progressEvent.loaded * 100) / totalLength)}%`)
        }
      },
      method: `post`
    })
      .then(response => {
        this._formData = new FormData()
      })
      .catch(error => {
        this._formData = new FormData()
        console.error(error)
      })
  }

  _cancelModal () {
    $(`#reference-file-name`).empty()
    $(`#bed-file-name`).empty()
    $(`#scheme-version`).val(``)
    $(`#scheme-name`).val(``)
    this._formData = new FormData()
  }

  _addEventsToDropZone () {
    const dropZone = $(`#drop-zone`)
    dropZone.on(`dragover`, event => {
      this._dragOverHandler(event)
    })
    dropZone.on(`dragleave`, event => {
      this._dragLeaveHandler(event)
    })
    dropZone.on(`dragenter`, event => {
      this._dragEnterHandler(event)
    })
    $(`#scheme-upload`).on(`change`, event => {
      console.log(`event trigg`)
      this._dropHandler(event, false)
    })
    $(`#cancel-scheme`).on(`click`, event => { this._cancelModal() })
  }

  _dragOverHandler (ev) {
    // Prevent default behavior (Prevent file from being opened)
    ev.preventDefault()
    $(ev.currentTarget).addClass(`drop-zone-entered`)
  }

  _dragEnterHandler (ev) {
    $(ev.currentTarget).addClass(`drop-zone-entered`)
  }

  _dragLeaveHandler (ev) {
    $(ev.currentTarget).removeClass(`drop-zone-entered`)
  }
}
