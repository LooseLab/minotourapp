class flowcellNotificationController {
    /**
     * Controller to handle the display of sent messages and the creation of conditions for messages
     *
     */
    constructor() {
        this._listMessages = [];
        this._messagesView = new MessagesView(document.querySelector("#all-messages"));
        this._messagesView.update(this._listMessages);
    }

}