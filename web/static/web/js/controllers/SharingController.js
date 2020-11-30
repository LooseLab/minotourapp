class SharingController {
  /**
   *
   * Create a sharing controller that is in charge of fetching and creating new user permissions on flowcells - used on the sharing-tab.html page
   */

  constructor (flowcellId) {
    this._flowcellId = flowcellId
    this._username = document.querySelector(`#id_username`)
    this._permission = document.querySelector(`#id_permission`)
    this._messageView = new MessageView(document.querySelector(`#sharing-messages`))
    this._dataTable = $(`#sharing_table`).DataTable({
      data: [],
      columns: [
        { data: `username` },
        { data: `permission` },
        { data: `user` }
      ],
      columnDefs: [
        {
          targets: 2,
          data: 2,
          render: (data, type, full, meta) => `<span href="#" onclick="mSharingController.delete('${full.permission_code}', ${full.user}, ${full.flowcell})">DELETE</span>`
        }]
    })
    this.getAll()
  }

  create (event) {
    event.preventDefault()
    const self = this
    const axiosInstance = axios.create({
      headers: { 'X-CSRFToken': getCookie(`csrftoken`) }
    })

    axiosInstance.post(`/api/v1/reads/lowcells/${this._flowcellId}/sharing/`, {
      flowcell: this._flowcellId,
      permission: this._permission.value,
      username: this._username.value
    }).then(response => {
      self.getAll()
      self._messageView.update(`New permission successfully created!`)
    }).catch(error => {
      self._messageView.update(`Something went wrong! Please check the following message: ${error.response.data.message}`)
    })
  }

  delete (permission, userId, flowcell) {
    const csrftoken = getCookie(`csrftoken`)
    const self = this
    const axiosInstance = axios.create({
      headers: { 'X-CSRFToken': csrftoken }
    })
    var url = `/api/v1/reads/flowcells/${this._flowcellId}/sharing/delete/`
    axiosInstance.post(url, {
      user_id: userId,
      permission
    }).then(response => {
      self.getAll()
      self._messageView.update(`Permission deleted!`)
    }).catch(error => {
      self._messageView.update(`Something went wrong! Please check the following message: ${error.response.data.message}`)
    })
  }

  getAll () {
    var self = this
    var url = `/api/v1/reads/flowcells/${this._flowcellId}/sharing/`

    axios.get(url)
      .then(response => {
        self._dataTable.clear().rows.add(response.data).draw()
      })
      .catch(error => {
        console.log(error)
      })
  }
}
