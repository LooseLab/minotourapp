class MTaskController {
  /* TODO rewrite this trash and redocument
    Create a task controller that is in charge of fetching and creating new tasks - used on the task-tab.html page,
     separate from Javascript that creates the drop-downs themselves
     */

  constructor (flowcell_id) {
    // store the flowcell id accessible to class instance
    this._flowcell_id = flowcell_id
    // The select box for the jobtype on the tasks tab page
    this._select_job_type = document.querySelector(`#id_job_type`)
    // The select box on the reference drop down on
    this._select_reference = document.querySelector(`#id_reference`)
    // target set value
    this._target_set = null
    // list for storing the tasks run or new tasks
    this._task_list = new MTaskList()
    // A instance of a message - getter and setter for its text
    this._message = new MMessage()
    // Create a MMessageView class for the message view element
    this._messageView = new MMessageView(document.querySelector(`#messageView`))
    // update the element with this._message
    this._messageView.update(this._message)
    // create a message view for reactivating the flowcell
    this._reactivatemessageView = new MMessageView(document.querySelector(`#messageViewReactivate`), `warning`)
    // update reactivate message
    // this._reactivatemessageView.update(this._message);
    // The name of the flowcell used for the deletion checks
    this._flowcellName = $(`#name-flowcell`).html()
    // Whether or not to delete the flowcell
    this._fromDatabase = $(`#from-database`)
  }

  _submitDeleteTask (taskNew) {
    // Not sure if this secure but ok
    const csrftoken = getCookie(`csrftoken`)
    const axiosInstance = axios.create({
      headers: { "X-CSRFToken": csrftoken }
    })

    axiosInstance.post(`/api/v1/tasks/`, {
      flowcell: taskNew.flowcell_id,
      reference: taskNew.reference_id,
      job_type: taskNew.job_type_id,
      target_set: taskNew.target_set_id
    }).then(function (response) {
      console.log(response)
      if (taskNew.job_type_id === `11`) {
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
  createNewJob (event) {
    let fromDatabase
    // prevent default behaviour
    event.preventDefault()
    this._reference_id = this._select_reference.value

    console.log(this._select_job_type.value)

    // If this is metagenomics - set the reference value to null and the target set id
    if (this._select_job_type.value === `10`) {
      // the target set text
      this._target_set = this._select_reference[this._select_reference.selectedIndex].text
      this._reference_id = null
    }

    // create an instance of the MTask class for this new task
    fromDatabase = this._fromDatabase[0].checked
    const task_new = new MTask(this._flowcell_id, this._select_job_type.value,
      this._reference_id, this._target_set, fromDatabase)
    // delete a flowcell, so show confirmation modal
    if (task_new.job_type_id === `11`) {
      $(`#confirm-submit`).modal()
      $(`#flowcell-name`).on(`input`, event => {
        const icon = $(`.fa-galactic-republic`)
        console.log(this._flowcellName)
        if (event.target.value === this._flowcellName) {
          console.log(`Great Success`)
          icon.removeClass(`disallowed`)
          icon.addClass(`allowed`)
          $(`#submit-taskdelete`).attr(`disabled`, false)
        } else {
          icon.removeClass(`allowed`)
          icon.addClass(`disallowed`)
          console.log(`Epic fail.`)
          $(`#submit-taskdelete`).attr(`disabled`, true)
        }
      })
      $(`#submit-taskdelete`).on(`click`, event => {
        $(`#confirm-submit`).modal(`hide`)
        this._submitDeleteTask(task_new)
      })
      return
    }

    // csrf safe
    const csrftoken = getCookie(`csrftoken`)
    // new request
    const xhr = new XMLHttpRequest()

    // Open the request for editing
    xhr.open(`POST`, `/api/v1/tasks/`)
    // set the request headers
    xhr.setRequestHeader(`X-CSRFToken`, csrftoken)
    xhr.setRequestHeader(`Content-Type`, `application/json`)
    xhr.setRequestHeader(`Accept`, `application/json`)

    const self = this

    xhr.onreadystatechange = function () {
      // event handler called in the readyState changing - when the request is submitted and returned
      if (this.readyState === XMLHttpRequest.DONE) {
        // if successful
        if (this.status === 200) {
          // If we have succesfully started a Delete flowcell task, return to the flowcell page
          // change the message to
          self._message.texto = `Task successfully created!`
          // update the message element to show the text
          self._messageView.update(self._message)
          // reload the task table to include the new task
          const taskTable = $(`.tasktable`)

          taskTable.DataTable().ajax.reload()
        } else {
          console.log(this)
          // something went wrong
          self._message.texto = `Something went wrong. Please check the following message. ` + this.responseText
          // update the message element to show the error text
          self._messageView.update(self._message)
        }
      }
    }
    // Send the xhr, with the data we need stringified
    xhr.send(JSON.stringify({

      flowcell: task_new.flowcell_id,
      reference: task_new.reference_id,
      job_type: task_new.job_type_id,
      target_set: task_new.target_set_id,
      from_database: task_new.fromDatabase
    }))
  }

  reactivateFlowcell (event) {
    event.preventDefault()

    const csrftoken = getCookie(`csrftoken`)
    // new request
    const xhr = new XMLHttpRequest()

    // Open the request for editing
    xhr.open(`POST`, `/api/v1/flowcells/reactivate/`)
    // set the request headers
    xhr.setRequestHeader(`X-CSRFToken`, csrftoken)
    xhr.setRequestHeader(`Content-Type`, `application/json`)
    xhr.setRequestHeader(`Accept`, `application/json`)
    const self = this
    xhr.onreadystatechange = function () {
      // event handler called in the readyState changing - when the request is submitted and returned
      if (this.readyState === XMLHttpRequest.DONE) {
        // if successful
        if (this.status === 200) {
          // change the message to
          self._message.texto = `Flowcell successfully reactivated`
          self._reactivatemessageView.update(self._message)
        } else {
          self._message.texto = this.responseText
          self._reactivatemessageView.update(self._message)
        }
      }
    }

    xhr.send(JSON.stringify({

      flowcell: this._flowcell_id

    }))
  }

  performActionOnTask (event, flowcellJobId, actionType) {
    const csrftoken = getCookie(`csrftoken`)
    // new request
    const xhr = new XMLHttpRequest()

    // Open the request for editing
    xhr.open(`POST`, `/api/v1/tasks/action`)

    // set the request headers
    xhr.setRequestHeader(`X-CSRFToken`, csrftoken)
    xhr.setRequestHeader(`Content-Type`, `application/json`)
    xhr.setRequestHeader(`Accept`, `application/json`)

    const self = this
    const lookup_action_type = { 1: `Reset`, 2: `Pause`, 3: `Delete` }
    const action = lookup_action_type[actionType]

    if (action === `Delete`) {
      const element = document.getElementById(`${action.toLowerCase()}_${flowcellJobId}`)
      element.style.visibility = `hidden`
    }

    self._message.texto = `${action}ing task with ID ${flowcellJobId} - please wait for confirmation`
    self._messageView.update(self._message)

    xhr.onreadystatechange = function () {
      // event handler called in the readyState changing - when the request is submitted and returned
      if (this.readyState === XMLHttpRequest.DONE) {
        // if successful
        if (this.status === 200) {
          // change the message to
          self._message.texto = this.responseText
          // update the message element to show the text
          self._messageView.update(self._message)

          // reload the task table to include the new task
          const taskTable = $(`.tasktable`)

          if (action === `Delete`) {
            console.log(`redraw later`)
          } else {
            taskTable.DataTable().ajax.reload()
          }
        } else {
          // something went wrong
          self._message.texto = `Something went wrong. Please check the following message. ` + this.responseText
          // update the message element to show the error text
          self._messageView.update(self._message)
        }
      }
    }
    // Send the xhr, with the data we need stringified
    xhr.send(JSON.stringify({
      flowcellJobId,
      actionType
    }))
  }
}
