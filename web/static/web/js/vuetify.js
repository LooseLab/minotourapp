// src/plugins/vuetify.js

// import {Vue} from '/static/web/libraries/vue.js';
// import Veutify from '/static/web/libraries/vuetify.js'

// import colors from 'web/node_modules/vuetify/lib/util/colors.js'

Vue.use(Veutify);

export default new Vuetify({
  theme: {
    themes: {
      light: {
        primary: "#E53935",
        secondary: "#FFCDD2",
        accent: "#3F51B5"
      },
    },
  },
})
