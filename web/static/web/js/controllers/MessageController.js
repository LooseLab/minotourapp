class MessageController {

    constructor() {
        this._messageTable = null;
    }

    drawMessageTable() {
        console.log("drawing table");
        $("#message-table").DataTable({
            ajax: {
                url: '/private/messages-sent',
                method: "GET",
            },
            columns: [
                {"data": "title"},
                {"data": "created_date"},
                {"data": "sender"},
                {"data": "content"},
            ]
        });
    }
}