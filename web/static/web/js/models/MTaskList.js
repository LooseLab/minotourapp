class MTaskList {
    /*
    Can contain either the current tasks ot the tasks type
     */
    constructor() {
        // list for storing the tasks
        this._tasks = [];
    }

    adiciona(task) {
        // Add a task to the list attached to the class
        this._tasks.push(task);
    }

    get tasks() {
        // return the tasks
        return [].concat(this._tasks);
    }
}