class MTaskList {

    constructor() {

        this._tasks = [];
    }

    adiciona(task) {

        this._tasks.push(task);
    }

    get tasks() {

        return [].concat(this._tasks);
    }
}