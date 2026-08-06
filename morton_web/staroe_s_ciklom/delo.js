

let funkaDlyaMortona = function(a,b){return 0;};
let wa=1;


function createMorton()
{
let xhr = new XMLHttpRequest();
xhr.open("GET", "morton.wasm");
xhr.responseType="arraybuffer";
xhr.onload =(e)=>
{
let qwe = e.target['response'];
let iO={"env":{}};
let wa = WebAssembly.instantiate(qwe,iO);
wa.then((e)=>{ wa = e["instance"]; funkaDlyaMortona=function(a,b){return wa.exports.morton(a,b); };});
};
xhr.send();

}




