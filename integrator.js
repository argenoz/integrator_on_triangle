
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
		e.set_points(ax,ay,bx,by,cx,cy);
		e.set_n(100,100);
		e.integrate_p();
		e.integrate_r();
		document.getElementById("resultat_p").innerText=e.get_Sp();
		document.getElementById("resultat_r").innerText=e.get_Sr();
		knopka=false;
	}

function delo()
{
	if(knopka)
		return;
	knopka=true;
	if(w==null)
		{
			let xhr = new XMLHttpRequest();
			xhr.open("GET","./integrator.wasm");
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
