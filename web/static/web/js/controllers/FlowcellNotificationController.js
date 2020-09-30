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
    this._referenceSet = new Set([`<option id="ref-place-holder" value="${-1}" selected>Please Choose</option>`])
    this._referenceSelect = $(`#reference`)
    this._contigSelect = $(`#chromosome`)
    this._table = $(`.extantNotif`)
    this.options = {
      backdrop: `true`
    }
    this._formData = new FormData()
    this._table = null
    this._disabledNotificationsHTML = {}
    this._type2Id = {
      arti: [`list-artic-top`, `list-artic`],
      suff: [`list-artic-top`, `list-artic`]
    }
    this._addListenerToForms()
  }

  /**
   * Function to call on reference select boxchange
   * @private
   */
  _referenceOnChange (event) {
    const referencePk = $(event.currentTarget).val()
    this._referenceSelect.find(`#ref-place-holder`).remove()
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
        const mappingOptions = $(`#list-mapping-top`)
        console.log(Object.keys(references))
        if (!Object.keys(references).length) {
          mappingOptions.addClass(`disabled`)
          mappingOptions.css(`cursor`, `not-allowed`)
          mappingOptions.attr(`title`, `No mapping tasks running`)
        } else {
          mappingOptions.removeClass(`disabled`)
          mappingOptions.css(`cursor`, `pointer`)
          mappingOptions.removeAttr(`title`)
        }
        Object.entries(references).forEach(([key, value]) => {
          console.log(key)
          console.log(value)
          // key is the reference name, value[0][2] is the reference PK
          this._referenceSet.add(`<option value="${value[0][2]}">${key}</option>`)
          // Value [0] is the contig name, value[1] is the contig PK
          value.forEach(([contigName, contigPk, referenceId]) => {
            if (!(referenceId in this._contigSet)) {
              this._contigSet[referenceId] = new Set()
              this._contigSet[referenceId].add(`<option id="place-holder" value="0" selected>Genome Coverage</option>`)
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
    $(`#multi-form-submit`).on(`click`, this.submitForms)
  }

  /**
   * Pull back extant notifications and remove them from the datalist, store their HTML in class, to be added if deleted
   * @private
   */
  _disableExtantNotifications () {
    const that = this
    this._axiosInstance.get(`/api/v1/communication/${this._flowcellId}/conditions/extant-types`).then(
      response => {
        const extantTypes = response.data
        if (extantTypes.includes(`arti`) && extantTypes.includes(`suff`)) {
          console.log(`hello`)
          if (!Object.prototype.hasOwnProperty.call(this._disabledNotificationsHTML, `list-artic-top`)) {
            that._disabledNotificationsHTML[`list-artic-top`] = $(`#list-artic-top`)[0].outerHTML
            $(`#list-artic-top`).remove()
            $(`#list-artic`).css(`display`, `none`)
          }
        }
        extantTypes.forEach(type => {
          console.log(type)
          console.log(![`arti`, `suff`].includes(type))
          const a = $(`#${type}-btn`)
          if (![`cov`, `arti`, `suff`].includes(type)) {
            if (!Object.prototype.hasOwnProperty.call(this._disabledNotificationsHTML, type)) {
              a.removeClass(`active`)
              that._disabledNotificationsHTML[type] = a[0].outerHTML
              a.remove()
            }
          } else if ([`arti`, `suff`].includes(type)) { $(`#${type}`).find(`#${type}`).prop(`disabled`, true) }
        })
        $(`#${$(`#list-tab-top`).children().first().attr(`id`)}`).tab(`show`)
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
    const formData = new FormData()
    event.preventDefault()

    console.log(event)
    $(event.currentTarget).find(`.noti-control`).each((index, field) => {
      console.log(index)
      console.log(field)
      const value = field.type === `checkbox` ? $(field).prop(`checked`) : $(field).val()
      if (value) {
        if (!formData.has(`flowcell`)) {
          formData.append(`notification_type`, $(event.currentTarget).attr(`id`))
          formData.append(`flowcell`, that._flowcellId)
        }
        formData.append(field.id.replace(`-`, `_`), value)
      }
    })
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
        that._disableExtantNotifications()
        that.getNotificationsForTable()
        $(event.currentTarget)[0].reset()
        that._formData = new FormData()
      }
    ).catch(
      error => {
        console.error(error)
        that._formData = new FormData()
      }
    )
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
          $(row).last().on(`click`, event => { this._deleteNotification(data.id, data.notification_type) })
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
    console.log(this._disabledNotificationsHTML)
    axiosInstance.delete(`/api/v1/communication/messages/conditions`, { params: { pk: notificationID } }).then(response => {
      const notificationType = response.data
      $(`.extantNotif`).DataTable().ajax.reload(null, false)
      if ([`arti`, `suff`].includes(notificationType)) {
        // append back thea artic tab
        if (!$(`#list-tab-tasks`).find(`#list-artic-top`).length) {
          $(`#list-tab-tasks`).append(this._disabledNotificationsHTML[`list-artic-top`])
          $(`#list-artic`).css(`display`, `block`)
          $(`#${notificationType}`).prop(`disabled`, false)
        }
      } else if (notificationType !== `cov`) {
        console.log(`WE ARE APP{EFNM`)
        $(`#list-tab-top`).append(this._disabledNotificationsHTML[notificationType])
      }
    }).catch(error => {
      console.error(error)
    })
  }

  // /**
  //  * Get the HTML
  //  * @private
  //  */
  // _getDataListHtml () {
  //
  // }

  /**
   * Render the page HTML anew. TODO this means we have to slightly change how we are doing this page
   */
  updateTab () {
    console.log(`updating tab`)
    this.getNotificationsForTable(this._flowcellId)
    this._disableExtantNotifications()
    this._getReferencesForCoverage()
  }
}
