class MSharingController {
    /*
    Create a sharing controller that is in charge of fetching and creating new user permissions on flowcells - used on the sharing-tab.html page
    */

    constructor(flowcell_id) {
        console.log('Initialising MSharingController');

        this._flowcell_id = flowcell_id;

        this._username = document.querySelector('#id_username');

        this._permission = document.querySelector('#id_permission');

        this._flowcelluserpermission_list = [];

        this._message = new MMessage();

        this._messageView = new MMessageView(document.querySelector('#sharing-messages'));

        this._messageView.update(this._message);

        this._dataTable = $('#sharing_table').DataTable({
            data: [],
            columns: [
                { data: 'username' },
                { data: 'permission' },
                { data: 'user' },
            ],
            columnDefs: [
                {
                    'targets': 2,
                    'data': 2,
                    'render': function(data, type, full, meta) {
                        return '<a href="#" onclick="mSharingController.delete(\'' + full['permission_code'] + '\', ' + full['user'] + ', ' + full['flowcell'] + ')">DELETE</a>';
                    }
                },]
        });

        this.getAll();
    }

    create(event) {
        event.preventDefault();

        let self = this;

        let csrftoken = getCookie('csrftoken');

        const axios_instance = axios.create({
            headers: {'X-CSRFToken': csrftoken}
        });

        axios_instance.post('/api/v1/flowcells/' + this._flowcell_id + '/sharing/', {
            flowcell: this._flowcell_id,
            permission: this._permission.value,
            username: this._username.value,
        }).then(function(response) {
            self.getAll();
            self._message.texto = 'New permission successfully created!';
            self._messageView.update(self._message);
        }).catch(function(error) {
            self._message.texto = 'Something went wrong! Please check the following message: ' + error.response.data.message;
            self._messageView.update(self._message);
        });
    }

    delete(permission, user_id, flowcell) {

        let csrftoken = getCookie('csrftoken');

        let self = this;

        const axios_instance = axios.create({
            headers: {'X-CSRFToken': csrftoken}
        });

        var url = '/api/v1/flowcells/' + this._flowcell_id + '/sharing/delete/';

        axios_instance.post(url, {
            user_id: user_id,
            permission: permission,
        }).then(function(response) {
            console.log(response);
            self.getAll();
            self._message.texto = 'Permission deleted!';
            self._messageView.update(self._message);
        }).catch(function(error) {
            self._message.texto = 'Something went wrong! Please check the following message: ' + error.response.data.message;
            self._messageView.update(self._message);
        });
    }

    getAll() {
        var self = this;
        var url = '/api/v1/flowcells/' + this._flowcell_id + '/sharing/';

        axios.get(url)
            .then(function(response) {
                console.log(self._dataTable.clear().rows.add(response.data).draw());
            })
            .catch(function(error) {
                console.log(error);
            });
    }
}
