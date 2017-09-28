var requestMessages = function () {

    var url_run = '/api/v1/messages';

    $.get(url_run, function (data) {
        console.log(data);
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

    });

    console.log('requesting new messages');
};

setInterval(function () {
    requestMessages();
}, 5000);


