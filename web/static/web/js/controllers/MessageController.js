class MessageController {

    constructor() {
        this._listMessages = [];
        this._messagesView = new MessagesView(document.querySelector("#all-messages"));
        this._messagesView.update(this._listMessages);
    }

    importMessages() {

        let service = new MessageService();
        service.getMessages((err, messages) => {

            if(err) {

                // this._mensagem.texto = err;
                return;
            }

            console.log(this);
            messages.forEach(message => this._listMessages.push(message));

            this._messagesView.update(this._listMessages);

            // this._mensagem.texto = 'Messages loaded with success.';
        });
    }

    requestMessages() {

        let duration = 30; // 30 seconds

    }

    template() {

        $('#all-messages').empty();

        for (var i = 0; i < data.length; i++) {
            $('#all-messages').append(
                $('<li>').append(
                    $('<a>').attr('href', '/messages/' + data[i].id).append(
                        $('<span>').attr('class', 'message').append(data[i].title)
                    )
                )
            );

            var messages_number = $('.message').length;
            $('#all-messages-number').html(messages_number);
        }
    }
}