class Flowcell {

    constructor(id=-1) {

        this._id = id;
        console.log(`Flowcell: ${this._id}`);
    }

    get id() {

        return this._id;
    }
}