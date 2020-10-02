class ReferenceManagementController {
  constructor () {
    this._axiosInstance = axios.create({
      headers: { 'X-CSRFToken': getCookie(`csrftoken`) }
    })
    // TODO make sure to refresh
    this._formData = new FormData()
    this._maxFileSize = 256000000
    this._datatableObj = $(`#reference-table`)
    this._addEventsToDropZone()
    this._createReferenceTable()
  }

  /**
   * Ajax request to delete a reference
   * @param referencePk {number} The primary key of the reference
   * @private
   */
  _deleteReference (referencePk) {
    this._axiosInstance.delete(`/api/v1/reference`, { data: { referencePk } }).then(
      response => {
        alert(response.data)
        this._datatableObj.DataTable().ajax.reload(null, false)
      }
    ).catch(
      error => {
        console.error(error.response)
        alert(`Reference failed to delete`)
      }
    )
  }

  /**
   * Create the reference Management dataTables table TODO update table function
   * @private
   */
  _createReferenceTable () {
    if ($.fn.DataTable.isDataTable(this._datatableObj)) {
      this._datatableObj.DataTable().ajax.reload(null, false)
    } else {
      this._datatableObj.DataTable({
        createdRow: (row, data, dataIndex, cells) => {
          $(row).find(`.ref-deletion`).on(`click`, (event) => {
            this._deleteReference(event.currentTarget.getAttribute(`data-delete`))
          })
        },
        ajax: {
          url: `/api/v1/reference`,
          async: true,
          error: (xhr, error, code) => {
            console.error(xhr)
            console.error(code)
          }
        },
        columnDefs: [
          {
            targets: 0,
            data: `name`

          },
          {
            targets: 1,
            data: `length`
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
    $(`#reference-upload`).on(`change`, event => {
      const formData = this._formData
      for (var i = 0; i < event.currentTarget.files.length; i++) {
        // If dropped items aren't files, reject them
        var file = event.currentTarget.files[i]
        if (file.size > this._maxFileSize) {
          this._modalError(`${file.name} is too large (${file.size})! Upload refused.`, 3000)
        }
        formData.append(`file_location`, file, file.name)
        formData.append(`file_name`, file.name)
      }
      // Post the files
      axios(
        {
          url: `/api/v1/reference/upload/`,
          data: formData,
          headers: {
            'X-CSRFToken': getCookie(`csrftoken`),
            'Content-Type': `multipart/form-data; boundary=${formData._boundary}`
          },
          onUploadProgress: progressEvent => {
            const totalLength = progressEvent.lengthComputable ? progressEvent.total : progressEvent.target.getResponseHeader(`content-length`) || progressEvent.target.getResponseHeader(`x-decompressed-content-length`)
            if (totalLength !== null) {
              $(`#uploadProgressBar`).css(`width`, `${Math.round((progressEvent.loaded * 100) / totalLength)}%`)
            }
            this._datatableObj.DataTable().ajax.reload(null, false)
          },
          method: `post`
        }
      )
        .then(
          response => {
            this._formData = new FormData()
            this._modalError(``, 0)
            this._datatableObj.DataTable().ajax.reload(null, false)
          })
        .catch(
          error => {
            console.error(error)
            if (error.status === 403) {
              this._modalError(error.response.data, 6000)
            }
            this._formData = new FormData()
          }
        )
    })
  }

  _dropHandler (ev) {
    // Prevent default behavior (Prevent file from being opened)
    ev.preventDefault()
    ev.stopImmediatePropagation()
    // check the md5 checksum here
    const formData = this._formData
    $(`#drop-zone`).removeClass(`drop-zone-entered`)
    $(`#mdbup`).html(`<div class="progress" id="uploadProgress">
          <div class="progress-bar progress-bar-striped progress-bar-animated" id="uploadProgressBar" role="progressbar" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100"></div>
        </div>`)
    if (ev.dataTransfer.items) {
      // Use DataTransferItemList interface to access the file(s)
      for (var i = 0; i < ev.dataTransfer.items.length; i++) {
        // If dropped items aren't files, reject them
        if (ev.dataTransfer.items[i].kind === `file`) {
          var file = ev.dataTransfer.items[i].getAsFile()
          if (file.size > this._maxFileSize) {
            this._modalError(`${file.name} is too large (${file.size})! Upload refused.`, 3000)
          }
          formData.append(`file_location`, file, file.name)
          formData.append(`file_name`, file.name)
        }
      }
      // Post the files
      axios(
        {
          url: `/api/v1/reference/upload/`,
          data: formData,
          headers: {
            'X-CSRFToken': getCookie(`csrftoken`),
            'Content-Type': `multipart/form-data; boundary=${formData._boundary}`
          },
          onUploadProgress: progressEvent => {
            const totalLength = progressEvent.lengthComputable ? progressEvent.total : progressEvent.target.getResponseHeader(`content-length`) || progressEvent.target.getResponseHeader(`x-decompressed-content-length`)
            if (totalLength !== null) {
              $(`#uploadProgressBar`).css(`width`, `${Math.round((progressEvent.loaded * 100) / totalLength)}%`)
            }
            this._datatableObj.DataTable().ajax.reload(null, false)
          },
          method: `post`
        }
      )
        .then(
          response => {
            this._formData = new FormData()
            this._modalError(``, 0)
            this._datatableObj.DataTable().ajax.reload(null, false)
          })
        .catch(
          error => {
            console.error(error)
            if (error.status === 403) {
              this._modalError(error.response.data, 6000)
            }
            this._formData = new FormData()
          }
        )
    } else {
      // Use DataTransfer interface to access the file(s)
      for (let i = 0; i < ev.dataTransfer.files.length; i++) {
        const file = ev.dataTransfer.items[i].getAsFile()
        formData.append(`file_location`, file, file.name)
        formData.append(`filename`, file.name)
      }
      axios(
        {
          url: `/api/v1/reference/upload/`,
          data: formData,
          headers: {
            'X-CSRFToken': getCookie(`csrftoken`),
            'Content-Type': `multipart/form-data; boundary=${formData._boundary}`
          },
          method: `post`
        }
      )
    }
  }

  /**
   *
   * @param message {string} Message to display
   * @param timeout {number} Time to wait until resetting modal
   * @private
   */
  _modalError (message, timeout) {
    if (message.length) {
      $(`#mdbup`).html(`<span style="font-size: 2rem">${message}</span>`)
    }
    setTimeout(() => {
      $(`#add-reference`).modal(`hide`)
      $(`#mdbup`).html(`<label id="drop-zone" class="dropzone" for="reference-upload">
                        Drag and drop or click to upload ...
                        </label>
                        <input type="file" id="reference-upload" style="display: none"
                           accept=".fsa,.fna,.fa,.fa.gz,.fsa.gz,.fna.gz" multiple>`)
      this._addEventsToDropZone()
    }, timeout)
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
