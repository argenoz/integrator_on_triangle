

let funkaDlyaMortona = function(a,b){return 0;};

let mort16_32 = function(a,b){return 0;};
let mort32_64 = function(a,b){return 0;};


function createMorton()
{
let xhr = new XMLHttpRequest();
xhr.open("GET", "morton_m.wasm");
xhr.responseType="arraybuffer";
xhr.onload =(e)=>
{
let qwe = e.target['response'];
let iO={"env":{}};
let wa = WebAssembly.instantiate(qwe,iO);
wa.then((e)=>{ 

wa = e["instance"]; mort16_32=function(a,b){

return wa.exports.morton16_32(a,b); };

mort32_64=function(a,b){

return wa.exports.morton32_64(a,b); };

}

);
};
xhr.send();

}




