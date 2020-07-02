class MessageController {
  // TODO add rewrite to interval update
  constructor () {
    this._messageTable = null
  }

  drawMessageTable () {
    console.log(`drawing table`)
    $(`#message-table`).DataTable({
      ajax: {
        url: `/web/private/sent_tweets/`,
        method: `GET`
      },
      columns: [
        { data: `title` },
        { data: `created_date` },
        { data: `sender_first_name` },
        { data: `flowcell_name` }
      ]
    })
  }
}
