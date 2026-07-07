
let flg=false;
let knopka=false;

let w=null;


let funka_=(x,y)=>{return 0;};

let funka=function(x,y){return funka_(x,y);};

function schet(x,y)
	{
		
		
	}


function vyzov()
	{
		let ax,ay,bx,by,cx,cy;
		let Su,tmp;
		let abso = (e)=>{if(e<0) e=-e;return e; };
		ax = document.getElementById('AX').value;
		ay = document.getElementById('AY').value;
		bx = document.getElementById('BX').value;
		by = document.getElementById('BY').value;
		cx = document.getElementById('CX').value;
		cy = document.getElementById('CY').value;
		funka_ = new Function('x','y',"{ return "+document.getElementById('funka').value+";}");
		//let v = w.exports.integrate(ax,ay,bx,by,cx,cy);
		//v = String(v);
		let e = w.exports;
		e.setV(ax,ay,bx,by,cx,cy);
		document.getElementById("resultat_p16").innerText = e.integrate_p16();
		document.getElementById("resultat_p3").innerText = e.integrate_p6();
		document.getElementById("resultat_p1").innerText = e.integrate_p3();
		knopka=false;
	}

function raschet()
{
	if(knopka)
		return;
	knopka=true;
	if(w==null)
		{
			let xhr = new XMLHttpRequest();
			xhr.open("GET","kod/integrator_all.wasm");
			xhr.responseType="arraybuffer";
			xhr.onload=(e)=>
				{
					let qwe = e.target["response"];
					let iO = {
						"env":{"f":funka}
						
					};
					let wa = WebAssembly.instantiate(qwe,iO);
					w = wa;
					w.then((e)=>{
								w = e["instance"];
								vyzov();
								});
					
					
					
				};
			xhr.send();
			return;
		}
		
	vyzov();	
	knopka=false;
}
let TTT = 0;

function showTheCode()
{
	let q = document.createElement("a");
	q.setAttribute("download","");
	q.setAttribute("href","kod/integrator_all.c");
	q.click();
	/*
	let q = document.createElement("div");
	if(document.getElementById("theCode")!=undefined) return;
	q.setAttribute("id","theCode");
	//document.body.appendChild(document.createElement("br"));
	//document.body.appendChild(q);
	document.getElementById("dlya_theCode").appendChild(q);
	let xhr = new XMLHttpRequest();
	xhr.open("GET","kod/integrator_all.c");
	//xhr.responseType="text";
	xhr.onload=(e)=>
	{
		qwe = document.createElement("code");
		qwe.innerText = e.target.response;
		document.getElementById("theCode").appendChild(qwe);
		
	};
	xhr.send();
	*/
}


