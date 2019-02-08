class DateHelper {

    constructor() {

        throw new Error('This class can not be instantiated');
    }

    static dateToText(date) {
        return `${date.getDate()}/${date.getMonth()+1}/${date.getFullYear()}`;
    }

    static datetimeToText(date) {
        return `${date.getDate()}/${date.getMonth()+1}/${date.getFullYear()} ${date.getHours()}:${date.getMinutes()}`;
    }

    static textToDate(text) {

        if(!/\d{4}-\d{2}-\d{2}/.test(text))
            throw new Error('Date must be in the format yyyy-mm-dd');

        return new Date(...text.split('-').map((item, indice) => item - indice % 2));
    }

}