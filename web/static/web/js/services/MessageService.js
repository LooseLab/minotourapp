class MessageService {

    getMessages(cb) {

        let xhr = new XMLHttpRequest();

        xhr.open('GET', '/api/v1/messages');

        xhr.onreadystatechange = function() {

            if(xhr.readyState == 4) {

                if(xhr.status == 200) {

                    console.log('Getting the messages.');
                    cb(null, JSON.parse(xhr.responseText)
                        .map(object => new Message(object.title, object.uuidstr, object.recipient_email, object.is_read, new Date(object.created_date))));
                } else {

                    console.log(xhr.responseText);
                    cb('It was not possible to get messages.', null);
                }
            }
        }

        xhr.send();
    }
}