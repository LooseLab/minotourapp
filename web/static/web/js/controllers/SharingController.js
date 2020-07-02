class SharingController {
  /* TODO rewrite and redocument all this
    Create a sharing controller that is in charge of fetching and creating new user permissions on flowcells - used on the sharing-tab.html page
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
          render: function (data, type, full, meta) {
            return `<a href="#" onclick="mSharingController.delete('` + full.permission_code + `', ` + full.user + `, ` + full.flowcell + `)">DELETE</a>`
          }
        }]
    })

    this.getAll()
  }

  create (event) {
    event.preventDefault()

    const self = this

    const csrftoken = getCookie(`csrftoken`)

    const axios_instance = axios.create({
      headers: { 'X-CSRFToken': csrftoken }
    })

    axios_instance.post(`/api/v1/flowcells/` + this._flowcellId + `/sharing/`, {
      flowcell: this._flowcell_id,
      permission: this._permission.value,
      username: this._username.value
    }).then(function (response) {
      self.getAll()
      self._messageView.update(`New permission successfully created!`)
    }).catch(function (error) {
      self._messageView.update(`Something went wrong! Please check the following message: ${error.response.data.message}`)
    })
  }

  delete (permission, userId, flowcell) {
    const csrftoken = getCookie(`csrftoken`)

    const self = this

    const axiosInstance = axios.create({
      headers: { 'X-CSRFToken': csrftoken }
    })

    var url = `/api/v1/flowcells/` + this._flowcellId + `/sharing/delete/`

    axiosInstance.post(url, {
      user_id: userId,
      permission
    }).then(function (response) {
      self.getAll()
      self._messageView.update(`Permission deleted!`)
    }).catch(function (error) {
      self._messageView.update(`Something went wrong! Please check the following message: ${error.response.data.message}`)
    })
  }

  getAll () {
    var self = this
    var url = `/api/v1/flowcells/` + this._flowcellId + `/sharing/`

    axios.get(url)
      .then(function (response) {
        self._dataTable.clear().rows.add(response.data).draw()
      })
      .catch(function (error) {
        console.log(error)
      })
  }
}
