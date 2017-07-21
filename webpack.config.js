var path = require('path');

module.exports = {
   entry: './web/static/web/minotour.js',
    output: {
        filename: 'bundle.js',
        path: path.resolve(__dirname, 'web/static/web'),
        library: 'Minotour',
        libraryTarget: 'var'
    }
};

