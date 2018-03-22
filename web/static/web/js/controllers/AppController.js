/*var doStuff = function () {

  	console.log('teste');
   	setTimeout(doStuff, appController._interval);
};

doStuff();*/

class AppController {

    constructor() {

        this._interval = 10000;
        this.updatePageContent();
    }

    updatePageContent() {

        console.log('ola');
        console.log(this);

    }
}