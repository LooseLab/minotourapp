class MessageController {
  // TODO add rewrite to interval update
  constructor () {
    this._messageTable = null
  }

  drawMessageTable () {
    console.log(`drawing table`)
    $(`#message-table`).DataTable({
      processing: true,
      serverSide: true,
      ajax: {
        url: `/api/v1/communication/sent-tweets/`,
        method: `GET`,
        async: true,
        error: (xhr, error, code) => {
          console.error(xhr)
          console.error(code)
        }
      },
      order: [[1, `desc`]],
      columns: [
        { data: `title` },
        { data: `created_date` },
        { data: `sender_first_name` },
        { data: `flowcell_name` }
      ]
    })
  }
}
