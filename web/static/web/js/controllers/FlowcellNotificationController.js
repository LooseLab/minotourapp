class FlowcellNotificationController {
  /**
   * Straight up just sorry about the code below
   * Summary. Controller to handle the display of sent messages and the creation of conditions for messages
   * @class flowcellID
   */
  constructor (flowcellId) {
    this._axiosInstance = axios.create({
      headers: { "X-CSRFToken": getCookie(`csrftoken`) }
    })
    this._flowcellId = flowcellId
    this._contigSet = {}
    this._referenceSet = new Set()
    this._referenceSelect = $(`#reference`)
    this._contigSelect = $(`#chromosome`)
    this._table = $(`.extantNotif`)
    this.options = {
      backdrop: `true`
    }
    this._formData = new FormData()
    this._table = null
    this.getNotificationsForTable(flowcellId)
    this._disableExtantNotifications()
    this._addListenerToForms()
    this._getReferencesForCoverage()
  }

  /**
   * Function to call on reference select boxchange
   * @private
   */
  _referenceOnChange (event) {
    const referencePk = $(event.currentTarget).val()
    console.log(event)
    console.log(referencePk)
    this._contigSelect.attr(`disabled`, false)
    this._contigSelect.html([...this._contigSet[referencePk]].join(``))
  }

  /**
   * Get the references and contigs for the coverage selects on condition create
   */
  _getReferencesForCoverage () {
    this._axiosInstance.get(`/api/v1/communication/${this._flowcellId}/conditions/references`).then(
      response => {
        const references = response.data
        console.log(references)
        Object.entries(references).forEach(([key, value]) => {
          console.log(key)
          console.log(value)
          // key is the reference name, value[0][2] is the reference PK
          this._referenceSet.add(`<option value="${value[0][2]}">${key}</option>`)
          // Value [0] is the contig name, value[1] is the contig PK
          value.forEach(([contigName, contigPk, referenceId]) => {
            console.log(contigPk)
            if (!(referenceId in this._contigSet)) {
              this._contigSet[referenceId] = new Set()
            }
            this._contigSet[referenceId].add(`<option value="${contigPk}">${contigName}</option>`)
          })
        })
        this._referenceSelect.html([...this._referenceSet].join(``))
      }
    )
  }

  _addListenerToForms () {
    $(`#reference`).on(`change`, event => { this._referenceOnChange(event) })
    $(`.notification-form`).on(`submit`, event => { this._postNotificationCreation(event, this) })
    $(`#multi-fomr-submit`).on(`click`, this.submitForms)
  }

  _disableExtantNotifications () {
    console.log(this)
    this._axiosInstance.get(`/api/v1/communication/${this._flowcellId}/conditions/extant-types`).then(
      response => {
        const extantTypes = response.data
        if (extantTypes.includes(`arti`) && extantTypes.includes(`suff`)) {
          console.log(`hello`)
          $(`#list-artic-top`).remove()
          $(`#list-artic`).remove()
        }
        extantTypes.forEach(type => {
          const a = $(`#${type}`)
          a.remove()
        })
        $(`#${$(`#list-tab`).children().first().attr(`id`)}`).tab(`show`)
      }
    ).catch(
      error => {
        console.error(error)
      }
    )
  }

  submitForms () {
    $(`#arti`).submit()
    $(`#suff`).submit()
  }

  _postNotificationCreation (event, that) {
    console.log($(event.currentTarget))
    console.log(that)
    const formData = new FormData()
    formData.append(`notification_type`, $(event.currentTarget).attr(`id`))
    formData.append(`flowcell`, that._flowcellId)
    console.log(formData)
    $(event.currentTarget).find(`.form-control`).each((index, field) => {
      console.log(field)
      that._formData.append(field.id.replace(`-`, `_`), $(field).val())
    })
    console.log(formData)
    axios({
      url: `/api/v1/communication/conditions/create-condition/`,
      data: formData,
      headers: {
        'X-CSRFToken': getCookie(`csrftoken`),
        'Content-Type': `multipart/form-data; boundary=${formData._boundary}`
      },
      method: `POST`
    }).then(
      response => {
        console.log(response)
        that._disableExtantNotifications()
        that.getNotificationsForTable()
        that._formData = new FormData()
      }
    ).catch(
      error => {
        console.error(error)
        that._formData = new FormData()
      }
    )

    console.log(event)
    event.preventDefault()
  }

  /**
   * @function Get all the Notifications we have already created for this flowcell.
   */
  getNotificationsForTable (flowcellId) {
    this._table = $(`.extantNotif`)
    if ($.fn.DataTable.isDataTable(this._table)) {
      this._table.DataTable().ajax.reload(null, false)
    } else {
      this._table.DataTable({
        createdRow: (row, data, index) => {
          $(row).addClass(`notification-row`)
          $(row).last().on(`click`, event => { this._deleteNotification(data.id) })
        },
        language: {
          infoEmpty: `No records available.`
        },
        ajax: {
          url: `/api/v1/communication/messages/conditions`,
          method: `GET`,
          data: { flowcellId }
        },
        columns: [
          { data: `id` },
          { data: `notification_type` },
          { data: `upper_limit` },
          { data: `lower_limit` },
          { data: `completed` },
          { data: `ref_name` },
          { data: `chrom_line_name` },
          { data: `coverage_target` },
          {
            data: null,
            orderable: false,
            width: `15%`,
            render: (data, type, full) => `<a id="delete_${data.id}"><i class="fa fa-times task-icon"></i> Delete </a>`
          }
        ]
      })
    }
  }

  /**
   * @function Delete a notification
   * @param notificationID - The PK of the notification to be deleted
   */
  _deleteNotification (notificationID) {
    const csrftoken = getCookie(`csrftoken`)

    const axiosInstance = axios.create({
      headers: { "X-CSRFToken": csrftoken }
    })
    console.log(notificationID)

    axiosInstance.delete(`/api/v1/communication/messages/conditions`, { params: { pk: notificationID } }).then(response => {
      console.log(response)
      $(`.extantNotif`).DataTable().ajax.reload(null, false)
    }).catch(error => {
      console.error(error)
    })
  }
}
