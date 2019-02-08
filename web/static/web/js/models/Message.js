class Message {

    constructor(title, uuidstr, recipient_email, is_read, created_date) {

        this._title = title;
        this._uuidstr = uuidstr;
        this._recipient_email = recipient_email;
        this._is_read = is_read;
        this._created_date = created_date;
    }

    get title() {

        return this._title;
    }

    get uuidstr() {

        return this._uuidstr;
    }

    get recipient_email() {

        return this._recipient_email;
    }

    get is_read() {

        return this._is_read;
    }

    get created_date() {

        return this._created_date;
    }

    /*set texto(texto) {

        this._texto = texto;
    }*/
}