/*! Built with http://stenciljs.com */
((w,d,x,n,h,c,r)=>{((s)=>{s&&(r=s.getAttribute('data-resources-url'))})(d.querySelector("script[data-namespace='docssite']"));
function e(e,t){return"sc-"+e.t+(t&&t!==l?"-"+t:"")}function t(e,t){return e+(t?"-h":"-s")}function o(e,t){let n,o,l=null,i=!1,s=!1,r=arguments.length;for(;r-- >2;)$.push(arguments[r]);for(;$.length>0;){let t=$.pop();if(t&&void 0!==t.pop)for(r=t.length;r--;)$.push(t[r]);else"boolean"==typeof t&&(t=null),(s="function"!=typeof e)&&(null==t?t="":"number"==typeof t?t=String(t):"string"!=typeof t&&(s=!1)),s&&i?l[l.length-1].vtext+=t:null===l?l=[s?{vtext:t}:t]:l.push(s?{vtext:t}:t),i=s}if(null!=t){if(t.className&&(t.class=t.className),"object"==typeof t.class){for(r in t.class)t.class[r]&&$.push(r);t.class=$.join(" "),$.length=0}null!=t.key&&(n=t.key),null!=t.name&&(o=t.name)}return"function"==typeof e?e(t,l||[],k):{vtag:e,vchildren:l,vtext:void 0,vattrs:t,vkey:n,vname:o,o:void 0,l:!1}}const l="$",i={},s={enter:13,escape:27,space:32,tab:9,left:37,up:38,right:39,down:40},a=(t,n,o,i)=>{let s=o.t+l,r=o[s];if((2===o.i||1===o.i&&!t.u.s)&&(i["s-sc"]=r?e(o,i.mode):e(o)),r){let e=n.p.head,o=t.m.get(e);if(o||t.m.set(e,o={}),!o[s]){let t;{t=r.content.cloneNode(!0),o[s]=!0;const l=e.querySelectorAll("[data-styles]");n.v(e,t,l.length&&l[l.length-1].nextSibling||e.firstChild)}}}},f=e=>null!=e,u=e=>e.toLowerCase(),p=(e,t,n,o,l,i)=>{if("class"!==n||i)if("style"===n){for(const e in o)l&&null!=l[e]||(/-/.test(e)?t.style.removeProperty(e):t.style[e]="");for(const e in l)o&&l[e]===o[e]||(/-/.test(e)?t.style.setProperty(e,l[e]):t.style[e]=l[e])}else if("o"!==n[0]||"n"!==n[1]||!/[A-Z]/.test(n[2])||n in t)if("list"!==n&&"type"!==n&&!i&&(n in t||-1!==["object","function"].indexOf(typeof l)&&null!==l)){const o=e.M(t);o&&o.g&&o.g[n]?m(t,n,l):"ref"!==n&&(m(t,n,null==l?"":l),null!=l&&!1!==l||e.u.k(t,n))}else null!=l&&"key"!==n?((e,t,n,o="boolean"==typeof n,l)=>{l=t!==(t=t.replace(/^xlink\:?/,"")),null==n||o&&(!n||"false"===n)?l?e.removeAttributeNS("http://www.w3.org/1999/xlink",u(t)):e.removeAttribute(t):"function"!=typeof n&&(n=o?"":n.toString(),l?e.setAttributeNS("http://www.w3.org/1999/xlink",u(t),n):e.setAttribute(t,n))})(t,n,l):(i||e.u.j(t,n)&&(null==l||!1===l))&&e.u.k(t,n);else n=u(n)in t?u(n.substring(2)):u(n[2])+n.substring(3),l?l!==o&&e.u.C(t,n,l,0):e.u.W(t,n,0);else if(o!==l){const e=b(o),n=b(l),i=e.filter(e=>!n.includes(e)),s=b(t.className).filter(e=>!i.includes(e)),r=n.filter(t=>!e.includes(t)&&!s.includes(t));s.push(...r),t.className=s.join(" ")}},b=e=>null==e||""===e?[]:e.trim().split(/\s+/),m=(e,t,n)=>{try{e[t]=n}catch(e){}},v=(e,t,n,o,l)=>{const s=11===n.o.nodeType&&n.o.host?n.o.host:n.o,r=t&&t.vattrs||i,a=n.vattrs||i;for(l in r)a&&null!=a[l]||null==r[l]||p(e,s,l,r[l],void 0,o,n.l);for(l in a)l in r&&a[l]===("value"===l||"checked"===l?s[l]:r[l])||p(e,s,l,r[l],a[l],o,n.l)};let y=!1;const M=(e,t)=>{e&&(e.vattrs&&e.vattrs.ref&&e.vattrs.ref(t?null:e.o),e.vchildren&&e.vchildren.forEach(e=>{M(e,t)}))},g=(e,t)=>{{let n=0,o=!1;const l=()=>t.performance.now(),i=!1!==e.asyncQueue,s=Promise.resolve(),r=[],a=[],c=[],f=[],u=t=>n=>{t.push(n),o||(o=!0,e.raf(b))},p=e=>{for(let t=0;t<e.length;t++)try{e[t](l())}catch(e){console.error(e)}e.length=0},d=(e,t)=>{let n,o=0;for(;o<e.length&&(n=l())<t;)try{e[o++](n)}catch(e){console.error(e)}o===e.length?e.length=0:0!==o&&e.splice(0,o)},b=()=>{n++,p(a);const t=i?l()+7*Math.ceil(n*(1/22)):Infinity;d(c,t),d(f,t),c.length>0&&(f.push(...c),c.length=0),(o=a.length+c.length+f.length>0)?e.raf(b):n=0};return e.raf||(e.raf=t.requestAnimationFrame.bind(t)),{tick(e){r.push(e),1===r.length&&s.then(()=>p(r))},read:u(a),write:u(c)}}},$=[],k={forEach:(e,t)=>e.forEach(t),map:(e,t)=>e.map(t)},j=(e,t,n)=>{const[o,l,,i,s,r]=e,a={color:{N:"color"}};if(i)for(t=0;t<i.length;t++)a[(n=i[t])[0]]={O:n[1],S:!!n[2],N:"string"==typeof n[3]?n[3]:n[3]?n[0]:0,A:n[4]};return{t:o,T:l,g:Object.assign({},a),i:s,R:r?r.map(C):void 0}},C=e=>({L:e[0],D:e[1],q:!!e[2],B:!!e[3],I:!!e[4]}),W=(e,t)=>f(t)&&"object"!=typeof t&&"function"!=typeof t?e===Boolean||4===e?"false"!==t&&(""===t||!!t):e===Number||8===e?parseFloat(t):e===String||2===e?t.toString():t:t,N=(e,t,n)=>{e.P.add(t),e.F.has(t)||(e.F.set(t,!0),e.H?e.queue.write(()=>O(e,t,n)):e.queue.tick(()=>O(e,t,n)))},O=async(e,n,l,i,s,r)=>{if(e.F.delete(n),!e.U.has(n)){if(!(s=e.Z.get(n))){if((r=e.G.get(n))&&!r["s-rn"])return void(r["s-rc"]=r["s-rc"]||[]).push(()=>{O(e,n,l)});if(s=L(e,n,e.J.get(n),l))try{s.componentWillLoad&&await s.componentWillLoad()}catch(t){e.K(t,3,n)}}((e,n,l,i)=>{try{const s=n.V.host,r=n.V.encapsulation,a=!1;let c,f=l;if(!l["s-rn"]){e.X(e,e.u,n,l);const o=l["s-sc"];o&&(e.u.Y(l,t(o,!0)),"scoped"===r&&e.u.Y(l,t(o)))}if(i.render||i.hostData||s||c){e._=!0;const t=i.render&&i.render();let n;e._=!1;const s=o(null,n,t),c=e.ee.get(l)||{};c.o=f,e.ee.set(l,e.render(l,c,s,a,r))}l["s-rn"]=!0,l["s-rc"]&&(l["s-rc"].forEach(e=>e()),l["s-rc"]=null)}catch(t){e._=!1,e.K(t,8,l,!0)}})(e,e.M(n),n,s),n["s-init"]()}},E=(e,t,n,o,l,i,s,r,a)=>{if(t.type||t.state){const c=e.te.get(n);t.state||(!t.attr||void 0!==c[l]&&""!==c[l]||(r=i&&i.ne)&&f(a=r[t.attr])&&(c[l]=W(t.type,a)),n.hasOwnProperty(l)&&(void 0===c[l]&&(c[l]=W(t.type,n[l])),"mode"!==l&&delete n[l])),o.hasOwnProperty(l)&&void 0===c[l]&&(c[l]=o[l]),t.watchCallbacks&&(c[R+l]=t.watchCallbacks.slice()),T(o,l,function c(t){return(t=e.te.get(e.oe.get(this)))&&t[l]},function u(n,o){(o=e.oe.get(this))&&(t.state||t.mutable)&&S(e,o,l,n,s)})}else if(t.elementRef)A(o,l,n);else if(t.context){const i=e.le(t.context);void 0!==i&&A(o,l,i.getContext&&i.getContext(n)||i)}},S=(e,t,n,o,l,i,s)=>{(s=e.te.get(t))||e.te.set(t,s={});const r=s[n];if(o!==r&&(s[n]=o,i=e.Z.get(t))){{const e=s[R+n];if(e)for(let t=0;t<e.length;t++)try{i[e[t]].call(i,o,r,n)}catch(e){console.error(e)}}!e._&&t["s-rn"]&&N(e,t,l)}},A=(e,t,n)=>{Object.defineProperty(e,t,{configurable:!0,value:n})},T=(e,t,n,o)=>{Object.defineProperty(e,t,{configurable:!0,get:n,set:o})},R="wc-",L=(e,t,n,o,l,i,s,r)=>{try{l=new(i=e.M(t).V),((e,t,n,o,l,i)=>{e.oe.set(o,n),e.te.has(n)||e.te.set(n,{}),Object.entries(Object.assign({color:{type:String}},t.properties,{mode:{type:String}})).forEach(([t,s])=>{E(e,s,n,o,t,l,i)})})(e,i,t,l,n,o),function a(e,t,n){if(t){const o=e.oe.get(n);t.forEach(t=>{n[t.method]={emit:n=>e.ie(o,t.name,{bubbles:t.bubbles,composed:t.composed,cancelable:t.cancelable,detail:n})}})}}(e,i.events,l);try{if(s=e.se.get(t)){for(r=0;r<s.length;r+=2)l[s[r]](s[r+1]);e.se.delete(t)}}catch(n){e.K(n,2,t)}}catch(n){l={},e.K(n,7,t,!0)}return e.Z.set(t,l),l},D=(e,t,n,o,l,i)=>{if(e.P.delete(t),(l=e.G.get(t))&&((o=l["s-ld"])&&((n=o.indexOf(t))>-1&&o.splice(n,1),o.length||l["s-init"]&&l["s-init"]()),e.G.delete(t)),e.re.length&&!e.P.size)for(;i=e.re.shift();)i()},q=(e,t,n,o)=>{t.forEach(([t,l])=>{3&l.O&&T(n,t,function n(){return(e.te.get(this)||{})[t]},function n(i){S(e,this,t,W(l.A,i),o)})})},B=(e,t,n,o,l)=>{if(n.connectedCallback=function(){((e,t,n)=>{e.ae.has(n)||(e.ae.set(n,!0),function o(e,t){const n=e.M(t);n.R&&n.R.forEach(n=>{n.q||e.u.C(t,n.L,function o(e,t,n,l){return o=>{(l=e.Z.get(t))?l[n](o):((l=e.se.get(t)||[]).push(n,o),e.se.set(t,l))}}(e,t,n.D),1,n.I,n.B)})}(e,n)),e.U.delete(n),e.ce.has(n)||(e.fe=!0,e.P.add(n),e.ce.set(n,!0),((e,t,n)=>{for(n=t;n=e.u.ue(n);)if(e.pe(n)){e.de.has(t)||(e.G.set(t,n),(n["s-ld"]=n["s-ld"]||[]).push(t));break}})(e,n),e.queue.tick(()=>{e.J.set(n,((e,t,n,o,l)=>(n.mode||(n.mode=e.be(n)),n["s-cr"]||e.me(n,"ssrv")||e.s&&1===t.i||(n["s-cr"]=e.ve(""),n["s-cr"]["s-cn"]=!0,e.v(n,n["s-cr"],e.ye(n)[0])),o={ne:{}},t.g&&Object.keys(t.g).forEach(i=>{(l=t.g[i].N)&&(o.ne[l]=e.me(n,l))}),o))(e.u,t,n)),e.he(t,n)}))})(e,t,this)},n.disconnectedCallback=function(){((e,t)=>{!e.we&&((e,t)=>{for(;t;){if(!e.Me(t))return 9!==e.ge(t);t=e.Me(t)}})(e.u,t)&&(e.U.set(t,!0),D(e,t),M(e.ee.get(t),!0),e.u.W(t),e.ae.delete(t),[e.G,e.$e,e.J].forEach(e=>e.delete(t)))})(e,this)},n["s-init"]=function(){((e,t,n,o,l,i,s)=>{if((l=e.Z.get(t))&&!e.U.has(t)&&(!t["s-ld"]||!t["s-ld"].length)){e.de.set(t,!0),(s=e.ke.has(t))||(e.ke.set(t,!0),t["s-ld"]=void 0,e.u.Y(t,n));try{M(e.ee.get(t)),(i=e.$e.get(t))&&(i.forEach(e=>e(t)),e.$e.delete(t)),!s&&l.componentDidLoad?l.componentDidLoad():s&&l.componentDidUpdate&&l.componentDidUpdate()}catch(n){e.K(n,4,t)}D(e,t)}})(e,this,o)},n.forceUpdate=function(){N(e,this,l)},t.g){const o=Object.entries(t.g);{let e={};o.forEach(([t,{N:n}])=>{n&&(e[n]=t)}),e=Object.assign({},e),n.attributeChangedCallback=function(t,n,o){(function l(e,t,n,o){const l=e[u(n)];l&&(t[l]=o)})(e,this,t,o)}}q(e,o,n,l)}};((e,t,n,i,r,c,p)=>{const d=n.performance,b={html:{}},m=n[e]=n[e]||{},h=((e,t,n)=>{const o=new WeakMap,l={p:n,s:!!n.documentElement.attachShadow,je:!1,ge:e=>e.nodeType,Ce:e=>n.createElement(e),We:(e,t)=>n.createElementNS(e,t),ve:e=>n.createTextNode(e),Ne:e=>n.createComment(e),v:(e,t,n)=>e.insertBefore(t,n),Oe:e=>e.remove(),xe:(e,t)=>e.appendChild(t),Y:(e,t)=>{e.classList.add(t)},ye:e=>e.childNodes,Me:e=>e.parentNode,Ee:e=>e.nextSibling,Se:e=>e.previousSibling,Ae:e=>u(e.nodeName),Te:e=>e.textContent,Re:(e,t)=>e.textContent=t,me:(e,t)=>e.getAttribute(t),Le:(e,t,n)=>e.setAttribute(t,n),k:(e,t)=>e.removeAttribute(t),j:(e,t)=>e.hasAttribute(t),be:t=>t.getAttribute("mode")||(e.Context||{}).mode,De:(e,o)=>"child"===o?e.firstElementChild:"parent"===o?l.ue(e):"body"===o?n.body:"document"===o?n:"window"===o?t:e,C:(t,n,i,r,a,c,f,u,p,d)=>{let b=t,m=i,v=o.get(t);d=n+r,v&&v[d]&&v[d](),"string"==typeof f?b=l.De(t,f):"object"==typeof f?b=f:(p=n.split(":")).length>1&&(b=l.De(t,p[0]),n=p[1]),b&&((p=n.split(".")).length>1&&(n=p[0],m=(e=>{e.keyCode===s[p[1]]&&i(e)})),u=l.je?{capture:!!a,passive:!!c}:!!a,e.ael(b,n,m,u),v||o.set(t,v={}),v[d]=(()=>{b&&e.rel(b,n,m,u),v[d]=null}))},W:(e,t,n,l)=>{(l=o.get(e))&&(t?l[t+n]&&l[t+n]():Object.keys(l).forEach(e=>{l[e]&&l[e]()}))},qe:(e,n,o,l)=>(l=new t.CustomEvent(n,o),e&&e.dispatchEvent(l),l),ue:(e,t)=>(t=l.Me(e))&&11===l.ge(t)?t.host:t,Be:(e,t,n,o)=>e.setAttributeNS(t,n,o)};e.ael||(e.ael=((e,t,n,o)=>e.addEventListener(t,n,o)),e.rel=((e,t,n,o)=>e.removeEventListener(t,n,o)));try{t.addEventListener("e",null,Object.defineProperty({},"passive",{get:()=>l.je=!0}))}catch(e){}return l})(m,n,i),w=h.p.documentElement,M=n["s-defined"]=n["s-defined"]||{},$=(e,t)=>{n.customElements.get(e.t)||(B(k,b[e.t]=e,t.prototype,c,d),t.observedAttributes=Object.values(e.g).map(e=>e.N).filter(e=>!!e),n.customElements.define(e.t,t))},k={u:h,Ie:$,M:e=>b[h.Ae(e)],le:e=>t[e],isClient:!0,pe:e=>!(!M[h.Ae(e)]&&!k.M(e)),K:(e,t,n)=>console.error(e,t,n&&n.tagName),queue:t.queue=g(m,n),he:(e,t)=>{{const n=e.T;let o=r+n+".entry.js";import(o).then(n=>{try{e.V=n[(e=>u(e).split("-").map(e=>e.charAt(0).toUpperCase()+e.slice(1)).join(""))(e.t)],function o(e,t,n,i,s){if(i){const n=t.t+(s||l);if(!t[n]){const o=e.Ce("template");t[n]=o,o.innerHTML=`<style>${i}</style>`,e.xe(e.p.head,o)}}}(h,e,e.i,e.V.style,e.V.styleMode),N(k,t,d)}catch(t){console.error(t),e.V=class{}}},e=>console.error(e,o))}},_:!1,H:!1,we:!1,X:a,G:new WeakMap,m:new WeakMap,ce:new WeakMap,ae:new WeakMap,ke:new WeakMap,de:new WeakMap,oe:new WeakMap,J:new WeakMap,Z:new WeakMap,U:new WeakMap,F:new WeakMap,$e:new WeakMap,se:new WeakMap,ee:new WeakMap,te:new WeakMap,P:new Set,re:[]};return t.isServer=t.isPrerender=!(t.isClient=!0),t.window=n,t.location=n.location,t.document=i,t.resourcesUrl=t.publicPath=r,t.enableListener=((e,t,n,o,l)=>(function i(e,t,n,o,l,s){if(t){const i=e.oe.get(t),r=e.M(i);if(r&&r.R)if(o){const o=r.R.find(e=>e.L===n);o&&e.u.C(i,n,e=>t[o.D](e),1,o.I,void 0===s?o.B:!!s,l)}else e.u.W(i,n,1)}})(k,e,t,n,o,l)),k.ie=t.emit=((e,n,o)=>h.qe(e,t.eventNameFn?t.eventNameFn(n):n,o)),m.h=o,m.Context=t,m.onReady=(()=>new Promise(e=>k.queue.write(()=>k.P.size?k.re.push(e):e()))),k.render=((e,t)=>{let n,o,l,i,s,r,a;const c=(l,p,d,b,m,h,w,M,g)=>{if(M=p.vchildren[d],n||(i=!0,"slot"===M.vtag&&(o&&t.Y(b,o+"-s"),M.vchildren?M.Pe=!0:M.Fe=!0)),f(M.vtext))M.o=t.ve(M.vtext);else if(M.Fe)M.o=t.ve("");else{if(h=M.o=y||"svg"===M.vtag?t.We("http://www.w3.org/2000/svg",M.vtag):t.Ce(M.Pe?"slot-fb":M.vtag),e.pe(h)&&e.de.delete(a),y="svg"===M.vtag||"foreignObject"!==M.vtag&&y,v(e,null,M,y),f(o)&&h["s-si"]!==o&&t.Y(h,h["s-si"]=o),M.vchildren)for(m=0;m<M.vchildren.length;++m)(w=c(l,M,m,h))&&t.xe(h,w);"svg"===M.vtag&&(y=!1)}return M.o["s-hn"]=r,(M.Pe||M.Fe)&&(M.o["s-sr"]=!0,M.o["s-cr"]=s,M.o["s-sn"]=M.vname||"",(g=l&&l.vchildren&&l.vchildren[d])&&g.vtag===M.vtag&&l.o&&u(l.o)),M.o},u=(n,o,l,s)=>{e.we=!0;const a=t.ye(n);for(l=a.length-1;l>=0;l--)(s=a[l])["s-hn"]!==r&&s["s-ol"]&&(t.Oe(s),t.v(h(s),s,m(s)),t.Oe(s["s-ol"]),s["s-ol"]=null,i=!0),o&&u(s,o);e.we=!1},p=(e,n,o,l,i,s,a,u)=>{const p=e["s-cr"];for((a=p&&t.Me(p)||e).shadowRoot&&t.Ae(a)===r&&(a=a.shadowRoot);i<=s;++i)l[i]&&(u=f(l[i].vtext)?t.ve(l[i].vtext):c(null,o,i,e))&&(l[i].o=u,t.v(a,u,m(n)))},d=(e,n,o,i)=>{for(;n<=o;++n)f(e[n])&&(i=e[n].o,l=!0,i["s-ol"]?t.Oe(i["s-ol"]):u(i,!0),t.Oe(i))},b=(e,t)=>e.vtag===t.vtag&&e.vkey===t.vkey&&("slot"!==e.vtag||e.vname===t.vname),m=e=>e&&e["s-ol"]?e["s-ol"]:e,h=e=>t.Me(e["s-ol"]?e["s-ol"]:e),w=(n,o,l)=>{const i=o.o=n.o,s=n.vchildren,r=o.vchildren;y=o.o&&f(t.ue(o.o))&&void 0!==o.o.ownerSVGElement,y="svg"===o.vtag||"foreignObject"!==o.vtag&&y,f(o.vtext)?(l=i["s-cr"])?t.Re(t.Me(l),o.vtext):n.vtext!==o.vtext&&t.Re(i,o.vtext):("slot"!==o.vtag&&v(e,n,o,y),f(s)&&f(r)?((e,n,o,l,i,s,r,a)=>{let v=0,y=0,M=n.length-1,g=n[0],$=n[M],k=l.length-1,j=l[0],C=l[k];for(;v<=M&&y<=k;)if(null==g)g=n[++v];else if(null==$)$=n[--M];else if(null==j)j=l[++y];else if(null==C)C=l[--k];else if(b(g,j))w(g,j),g=n[++v],j=l[++y];else if(b($,C))w($,C),$=n[--M],C=l[--k];else if(b(g,C))"slot"!==g.vtag&&"slot"!==C.vtag||u(t.Me(g.o)),w(g,C),t.v(e,g.o,t.Ee($.o)),g=n[++v],C=l[--k];else if(b($,j))"slot"!==g.vtag&&"slot"!==C.vtag||u(t.Me($.o)),w($,j),t.v(e,$.o,g.o),$=n[--M],j=l[++y];else{for(i=null,s=v;s<=M;++s)if(n[s]&&f(n[s].vkey)&&n[s].vkey===j.vkey){i=s;break}f(i)?((a=n[i]).vtag!==j.vtag?r=c(n&&n[y],o,i,e):(w(a,j),n[i]=void 0,r=a.o),j=l[++y]):(r=c(n&&n[y],o,y,e),j=l[++y]),r&&t.v(h(g.o),r,m(g.o))}v>M?p(e,null==l[k+1]?null:l[k+1].o,o,l,y,k):y>k&&d(n,v,M)})(i,s,o,r):f(r)?(f(n.vtext)&&t.Re(i,""),p(i,null,o,r,0,r.length-1)):f(s)&&d(s,0,s.length-1)),y&&"svg"===o.vtag&&(y=!1)},M=(e,n,o,l,i,s,r,a)=>{for(l=0,i=(o=t.ye(e)).length;l<i;l++)if(n=o[l],1===t.ge(n)){if(n["s-sr"])for(r=n["s-sn"],n.hidden=!1,s=0;s<i;s++)if(o[s]["s-hn"]!==n["s-hn"])if(a=t.ge(o[s]),""!==r){if(1===a&&r===t.me(o[s],"slot")){n.hidden=!0;break}}else if(1===a||3===a&&""!==t.Te(o[s]).trim()){n.hidden=!0;break}M(n)}},g=[],$=(e,n,o,i,s,r,a,c,f,u)=>{for(s=0,r=(n=t.ye(e)).length;s<r;s++){if((o=n[s])["s-sr"]&&(i=o["s-cr"]))for(c=t.ye(t.Me(i)),f=o["s-sn"],a=c.length-1;a>=0;a--)(i=c[a])["s-cn"]||i["s-nr"]||i["s-hn"]===o["s-hn"]||((3===(u=t.ge(i))||8===u)&&""===f||1===u&&null===t.me(i,"slot")&&""===f||1===u&&t.me(i,"slot")===f)&&(g.some(e=>e.He===i)||(l=!0,i["s-sn"]=f,g.push({Qe:o,He:i})));1===t.ge(o)&&$(o)}};return(c,f,u,p,d,b,m,v,y,h,k,j)=>{if(a=c,r=t.Ae(a),s=a["s-cr"],n=p,o=a["s-sc"],i=l=!1,w(f,u),i){for($(u.o),m=0;m<g.length;m++)(v=g[m]).He["s-ol"]||((y=t.ve(""))["s-nr"]=v.He,t.v(t.Me(v.He),v.He["s-ol"]=y,v.He));for(e.we=!0,m=0;m<g.length;m++){for(v=g[m],k=t.Me(v.Qe),j=t.Ee(v.Qe),y=v.He["s-ol"];y=t.Se(y);)if((h=y["s-nr"])&&h&&h["s-sn"]===v.He["s-sn"]&&k===t.Me(h)&&(h=t.Ee(h))&&h&&!h["s-nr"]){j=h;break}(!j&&k!==t.Me(v.He)||t.Ee(v.He)!==j)&&v.He!==j&&(t.Oe(v.He),t.v(k,v.He,j))}e.we=!1}return l&&M(u.o),g.length=0,u}})(k,h),w["s-ld"]=[],w["s-rn"]=!0,w["s-init"]=(()=>{k.de.set(w,m.loaded=k.H=!0),h.qe(n,"appload",{detail:{namespace:e}})}),p.map(j).forEach(e=>$(e,class extends HTMLElement{})),k.fe||w["s-init"](),((e,t,n,o,l,i)=>{if(t.componentOnReady=((t,n)=>{if(!t.nodeName.includes("-"))return n(null),!1;const o=e.M(t);if(o)if(e.de.has(t))n(t);else{const o=e.$e.get(t)||[];o.push(n),e.$e.set(t,o)}return!!o}),l){for(i=l.length-1;i>=0;i--)t.componentOnReady(l[i][0],l[i][1])&&l.splice(i,1);for(i=0;i<o.length;i++)if(!n[o[i]].componentOnReady)return;for(i=0;i<l.length;i++)l[i][1](null);l.length=0}})(k,m,n,n["s-apps"],n["s-cr"]),m.initialized=!0,k})(n,x,w,d,r,h,c);
})(window,document,{},"DocsSite","hydrated",[["footer-bar","mrv57at4",1],["header-bar","ngcitrgw",1,[["el",64],["isMobileMenuShown",16],["isSearchVisible",1,0,"is-search-visible",4],["isSticky",16],["query",1,0,1,2],["version",1,0,1,2]],0,[["window:scroll","handleScroll",0,1],["window:resize","handleResize",0,1]]],["icon-list","zdmwt86j",1,[["data",1,0,1,1],["el",64],["isHeaderSearchVisible",16],["query",1,0,1,2],["selectedIcon",16],["selectedIconType",16]],0,[["body:keyup","escListener"],["body:click","handleBodyClicked"],["clearToast","handleClearToast"],["window:scroll","handleScroll",0,1]]],["icon-search","pnozjohf",1,[["autofocus",1,0,1,2],["query",1,0,1,2],["showClearCtrl",16],["size",1,0,1,2]],0,[["keyup","searchListener"]]],["ionicons-site","ngcitrgw",1,[["data",16],["isHeaderSearchVisible",16],["query",16]],0,[["window:scroll","handleScroll",0,1],["hasSearched","searchHandler"],["toggleHeaderSearch","toggleHandler"]]],["landing-page","zdmwt86j",1,[["data",1,0,1,1],["el",64],["query",1,0,1,2]]],["notfound-page","8j0j37e6",1],["stencil-route","ngcitrgw",1,[["component",1,0,1,2],["componentProps",1],["componentUpdated",1],["el",64],["exact",1,0,1,4],["group",1,0,1,2],["history",1],["historyType",1,0,"history-type",2],["location",1],["match",2],["routeRender",1],["routeViewsUpdated",1],["scrollTopOffset",1,0,"scroll-top-offset",8],["url",1,0,1,2]]],["stencil-route-link","m3zdr5uv",0,[["activeClass",1,0,"active-class",2],["anchorClass",1,0,"anchor-class",2],["anchorId",1,0,"anchor-id",2],["anchorRole",1,0,"anchor-role",2],["anchorTabIndex",1,0,"anchor-tab-index",2],["anchorTitle",1,0,"anchor-title",2],["ariaHaspopup",1,0,"aria-haspopup",2],["ariaLabel",1,0,"aria-label",2],["ariaPosinset",1,0,"aria-posinset",2],["ariaSetsize",1,0,"aria-setsize",8],["custom",1,0,1,2],["el",64],["exact",1,0,1,4],["history",1],["location",1],["match",16],["root",1,0,1,2],["strict",1,0,1,4],["url",1,0,1,2],["urlMatch",1,0,"url-match",2]]],["stencil-route-switch","ngcitrgw",0,[["el",64],["group",1,0,1,2],["location",1],["queue",4,0,0,0,"queue"],["routeViewsUpdated",1],["scrollTopOffset",1,0,"scroll-top-offset",8]]],["stencil-router","ngcitrgw",0,[["history",16],["historyType",1,0,"history-type",2],["isServer",4,0,0,0,"isServer"],["location",16],["queue",4,0,0,0,"queue"],["root",1,0,1,2],["scrollTopOffset",1,0,"scroll-top-offset",8],["titleSuffix",1,0,"title-suffix",2]]],["toast-bar","zdmwt86j",1,[["el",64],["hadIconOnce",16],["selectedIcon",1],["selectedIconType",1,0,"selected-icon-type",2],["showCopiedConfirm",16],["touchEndY",16],["touchStartY",16]]],["usage-page","jkivsoqg",1,[["data",1,0,1,1],["exampleIcon",16],["exampleType",16],["match",1],["queue",4,0,0,0,"queue"]]]]);