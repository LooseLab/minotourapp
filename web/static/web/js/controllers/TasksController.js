class TasksController {
  /* TODO rewrite this trash and redocument
    Create a task controller that is in charge of fetching and creating new tasks - used on the task-tab.html page,
     separate from Javascript that creates the drop-downs themselves
     */

  constructor (flowcellId) {
    // store the flowcell id accessible to class instance
    this._flowcellId = flowcellId
    // The select box for the jobtype on the tasks tab page
    this._selectJobType = $(`#job-type-select`)
    // The select box on the reference drop down on
    this._selectReference = $(`#reference-select`)
    // target set value
    this._target_set = null
    // Create a MessageView class for the message view element
    this._messageView = new MessageView(document.querySelector(`#messageView`))
    // create a message view for reactivating the flowcell
    this._reactivatemessageView = new MessageView(document.querySelector(`#messageViewReactivate`), `warning`)
    // update reactivate message
    // The name of the flowcell used for the deletion checks
    this._flowcellName = $(`#name-flowcell`).html()
    // Whether or not to delete the flowcell
    this._fromDatabase = $(`#from-database`)
    this._axiosInstance = axios.create({
      headers: { "X-CSRFToken": getCookie(`csrftoken`) }
    })
    this._selectOptionsList = new Set()
    this._selectOptionsList.add(`<option value="-1">-- Please Choose --</option>`)
    this._referenceSelectOptionsList = new Set()
    this._referenceSelectOptionsList.add(`<option value="-1">-- Please Choose --</option>`)
    this._targetSetList = new Set()
    this._addListenerToTasksForm()
    this._addListenerToReactivateButton()
    this._targetSetList.add(`<option value="-1">-- Please Choose --</option>`)
    this._interval = setInterval(this._flowcellTaskHistoryTable, 30000, this._flowcellId)
    $(window).on(`unload`, function () {
      console.log(`clearing task interval`)
      clearInterval(this._interval)
    })
  }

  _addListenerToReactivateButton () {
    $(`#reactivate_flowcell`).on(`click`, () => {
      this.reactivateFlowcell(event)
    })
  }

  _addListenerToTasksForm () {
    $(`#post-form-task-create`).on(`submit`, () => {
      this._createNewJob(event)
    })
  }

  updateTab () {
    this._flowcellTaskHistoryTable(this._flowcellId)
    this._loadTasksForm()
  }

  /**
   * Submit a delete flowcell task specifically. Written seperately as it requires additional checks first before
   * submission
   * @param taskNew {MinotourTask} A Minotourtask instance
   * @private
   */
  _submitDeleteTask (taskNew) {
    this._axiosInstance.post(`/api/v1/tasks/`, {
      flowcell: taskNew.flowcellId,
      reference: taskNew.referenceId,
      job_type: taskNew.jobTypeId,
      target_set: taskNew.targetSetId
    }).then(function (response) {
      console.log(response)
      if (taskNew.jobTypeId === `11`) {
        window.location.href = `/`
        $(`#basicModal`).modal()
        setTimeout(function () {
          $(`#basicModal`).modal(`hide`)
        }, 6000)
      }
    }).catch(function (response) {
      $(`#basicModal`).modal()
    })
  }

  // Create a new event handler for when we submit a new job
  _createNewJob (event) {
    const fromDatabase = this._fromDatabase[0].checked
    const referenceId = this._selectReference.val() === `-1` ? null : this._selectReference.val()
    const jobTypeId = this._selectJobType.val()
    const targetSetId = this._selectReference.val()
    // create an instance of the MinotourTask class for this new task
    const taskNew = new MinotourTask(this._flowcellId, jobTypeId, referenceId,
      targetSetId, fromDatabase)
    const taskTable = $(`.tasktable`)
    const self = this

    // prevent default behaviour
    event.preventDefault()
    // delete a flowcell, so show confirmation modal
    if (taskNew.jobTypeId === `11`) {
      $(`#confirm-submit`).modal(`show`)
      $(`#flowcell-name`).on(`input`, event => {
        const icon = $(`.fa-galactic-republic`)
        console.log(this._flowcellName)
        if (event.target.value === this._flowcellName) {
          icon.removeClass(`disallowed`)
          icon.addClass(`allowed`)
          $(`#submit-task-delete`).attr(`disabled`, false)
        } else {
          icon.removeClass(`allowed`)
          icon.addClass(`disallowed`)
          $(`#submit-task-delete`).attr(`disabled`, true)
        }
      })
      $(`#submit-task-delete`).on(`click`, event => {
        $(`#confirm-submit`).modal(`hide`)
        this._submitDeleteTask(taskNew)
      })
      return
    }

    this._axiosInstance.post(`/api/v1/tasks/`, {
      flowcell: taskNew.flowcellId,
      reference: taskNew.referenceId,
      job_type: taskNew.jobTypeId,
      target_set: taskNew.targetSetId,
      from_database: taskNew.fromDatabase
    })
      .then(
        response => {
          self._messageView.update(`Task successfully created!`)
          // reload the task table to include the new task
          taskTable.DataTable().ajax.reload(null, false)
        }
      ).catch(
        error => {
          console.error(error)
          self._messageView.update(`Problem creating task: ${error.response.data.error_messages.non_field_errors[0]}`)
        }
      )
  }

  /**
   * Reactivate a flowcell record in the database.
   * @param event
   */
  reactivateFlowcell (event) {
    event.preventDefault()
    this._axiosInstance.post(`/api/v1/flowcells/${this._flowcellId}/reactivate/`).then(
      response => {
        this._reactivatemessageView.update(`Flowcell successfully reactivated`)
      }
    ).catch(
      error => {
        console.error(error)
        this._reactivatemessageView.update(`Issue updating flowcell - ${error.response.data}`)
      }
    )
  }

  performActionOnTask (event, flowcellJobId, actionType) {
    const self = this
    const lookupActionType = { 1: `Reset`, 2: `Pause`, 3: `Delete` }
    const action = lookupActionType[actionType]
    const taskTable = $(`.tasktable`)

    if (action === `Delete`) {
      document.getElementById(`${action.toLowerCase()}_${flowcellJobId}`).style.visibility = `hidden`
    }
    self._messageView.update(`${action}ing task with ID ${flowcellJobId} - please wait for confirmation`)
    this._axiosInstance.post(
      `/api/v1/tasks/action`,
      {
        flowcellJobId,
        actionType
      }
    ).then(
      response => {
        self._messageView.update(response.data)
        taskTable.DataTable().ajax.reload(null, false)
      }
    ).catch(
      error => {
        console.error(error)
        self._messageView.update(`Something went wrong. Please check the following message: ${error.response.data}`)
      }
    )
  }

  _flowcellTaskHistoryTable (flowcellId) {
  // Draw the Task table on the tasks tab showing previous and current analyses
    const table = $(`.tasktable`)
    // if it's already loaded don't reiniiialise
    if ($.fn.DataTable.isDataTable(table)) {
      table.DataTable().ajax.reload(null, false)
    } else {
      table.DataTable({
        ajax: {
          url: `/api/v1/tasks/?search_criteria=flowcell&search_value=${flowcellId.toString()}`,
          method: `GET`,
          async: true,
          error: (xhr, error, code) => {
            console.error(xhr)
            console.error(code)
          }
        },
        columns: [
          { data: `id` },
          { data: `task_type_name` },
          { data: `read_count` },
          { data: `running` },
          { data: `complete` },
          { data: `reference_name` },
          { data: `server_initiated` },
          { data: `from_database` },
          {
            data: null,
            orderable: false,
            width: `15%`,
            render: (data, type, full) => {
              if (!data.server_initiated) {
                return [`<a class="btn icon-task"  id="pause_${data.id}" onclick="flowcellController.flowcellTabController.tasksController.performActionOnTask(event, ${data.id}, 2)" ><i class="fa fa-${data.icon} task-icon"></i> ${data.iconText} </a>·
                                <a class="btn icon-task" id="reset_${data.id}" onclick="flowcellController.flowcellTabController.tasksController.performActionOnTask(event, ${data.id}, 1)"><i class="fa fa-recycle task-icon"></i> Reset </a>·
                                <a class="btn icon-task" id="delete_${data.id}" onclick="flowcellController.flowcellTabController.tasksController.performActionOnTask(event, ${data.id}, 3)"><i class="fa fa-times task-icon"></i> Delete </a>`]
              } else {
                return [`<a></a>`]
              }
            }
          }
        ]
      })
    }
  }

  _loadTasksForm () {
  // Select the div for the dropdown
    const referenceSelect = this._selectReference
    // add an event listener for change in the dropdown
    this._selectJobType.on(`change`, event => {
      const jobTypeId = this._selectJobType.val()
      // If there is not a value set
      const label = $(`#select-label`)
      if (jobTypeId === `10`) {
        label.html(`Target Sets`)
        axios.get(`/api/v1/metagenomics/targetsets`).then(
          response => {
            response.data.forEach(targetSet => {
              this._targetSetList.add(`<option value="${targetSet.target_set}">${targetSet.target_set}</option>`)
            })
            referenceSelect.html([...this._targetSetSet].join(``))
          }
        )
        referenceSelect.attr(`disabled`, false)
      } else if ([`4`, `16`].includes(jobTypeId)) {
        referenceSelect.attr(`disabled`, false)
      } else {
        referenceSelect.attr(`disabled`, true)
        if (label.html() !== `Reference`) {
          label.html(`Reference`)
          axios.get(`/api/v1/reference/`).then(
            response => {
              const references = response.data.data
              references.forEach(reference => {
                this._referenceSelectOptionsList.add(`<option value="${reference.id}">${reference.name}</option>`)
              })
              referenceSelect.html([...this._referenceSelectOptionsList].join(``))
            }
          )
        }
      }
    })
    // url to fetch all the available task types
    // fetch all available tasks types
    axios.get(`/api/v1/tasktypes/`).then(
      response => {
        // get the data
        const data = response.data.data
        // create option placeholder
        data.forEach(jobType => {
          this._selectOptionsList.add(`<option value="${jobType.id}">${jobType.name}</option>`)
        })
        this._selectJobType.html([...this._selectOptionsList].join(``))
      })
    // Do the same with the reference
    axios.get(`/api/v1/reference/`).then(
      response => {
        const references = response.data.data
        references.forEach(reference => {
          this._referenceSelectOptionsList.add(`<option value="${reference.id}">${reference.name}</option>`)
        })
        referenceSelect.html([...this._referenceSelectOptionsList].join(``))
      }
    )
  }
}
