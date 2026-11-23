/* Ideal-gas quasi-1D nozzle demonstration. Independent of the eductor closure. */
const Compressible=(()=>{
 const gamma=1.4,R=287.05,cp=gamma*R/(gamma-1);
 const pr=M=>(1+(gamma-1)*M*M/2)**(-gamma/(gamma-1));
 const area=M=>(2/(gamma+1)*(1+(gamma-1)*M*M/2))**((gamma+1)/(2*(gamma-1)))/M;
 function bisect(fn,lo,hi){let fl=fn(lo);for(let i=0;i<85;i++){let m=(lo+hi)/2,f=fn(m);if(fl*f<=0)hi=m;else{lo=m;fl=f;}}return (lo+hi)/2;}
 function mach(a,branch='sub'){if(a<1-1e-9)throw Error('Area below sonic area.');if(a<=1+1e-10)return 1;return branch==='sub'?bisect(m=>area(m)-a,1e-8,1):bisect(m=>area(m)-a,1,40);}
 function shock(M){if(M<1)throw Error('A normal shock requires upstream Mach ≥ 1.');const M2=Math.sqrt((1+(gamma-1)*M*M/2)/(gamma*M*M-(gamma-1)/2)),pRatio=1+2*gamma/(gamma+1)*(M*M-1),rhoRatio=(gamma+1)*M*M/(2+(gamma-1)*M*M),p0Ratio=pRatio*pr(M)/pr(M2);return {M2,pRatio,rhoRatio,tRatio:pRatio/rhoRatio,p0Ratio};}
 function state(M,p0,T0){const T=T0/(1+(gamma-1)*M*M/2),p=p0*pr(M),rho=p/(R*T),u=M*Math.sqrt(gamma*R*T);return {M,p,T,rho,u,p0,T0};}
 function assess({p0,pb,ratio}){
  if(![p0,pb,ratio].every(Number.isFinite)||p0<=0||pb<=0||ratio<1||ratio>20)return {needed:false,valid:false,reason:'Enter positive pressures and an area ratio from 1 to 20.'};
  if(pb>=p0)return {needed:false,valid:false,reason:'No forward pressure drive; no nozzle calculation is run.'};
  const jetMach=Math.sqrt(2/(gamma-1)*((p0/pb)**((gamma-1)/gamma)-1)),sonicRatio=pr(mach(ratio)),choked=pb/p0<=sonicRatio;
  const throatMach=choked?1:mach(area(jetMach)/ratio),peakMach=Math.max(jetMach,throatMach),needed=peakMach>=.3;
  return {needed,valid:true,choked,jetMach,throatMach,peakMach,reason:choked?'Activated: the ideal nozzle throat is choked.':needed?'Activated: estimated peak nozzle Mach is at least 0.30.':'Inactive: estimated peak nozzle Mach is below 0.30; retain the low-Mach treatment.'};
 }
 function solve({p0,pb,T0,diameter,ratio}){
  if(![p0,pb,T0,diameter,ratio].every(Number.isFinite)||p0<=0||pb<=0||T0<=0||diameter<=0||ratio<1||ratio>20)throw Error('Positive pressures, temperature, diameter and area ratio 1–20 are required.');
  if(pb>=p0)return {ok:false,reason:'No forward nozzle solution: backpressure must be below reservoir total pressure.'};
  const At=Math.PI*diameter**2/4,ms=mach(ratio),mu=mach(ratio,'sup'),onset=pr(ms),design=pr(mu),exitShock=design*shock(mu).pRatio,b=pb/p0;
  let Astar=At,shockArea=null,p02=p0,Astar2=null,regime,choked=b<=onset,Me;
  if(!choked){Me=Math.sqrt(2/(gamma-1)*((p0/pb)**((gamma-1)/gamma)-1));Astar=At*ratio/area(Me);regime='Unchoked · subsonic throughout';}
  else if(ratio<=1+1e-8){Me=1;regime='Choked · converging nozzle';}
  else if(b>exitShock){
   const exitAt=a=>{const s=shock(mach(a,'sup')),a2=a/area(s.M2),me=mach(ratio/a2);return s.p0Ratio*pr(me);};
   shockArea=bisect(a=>exitAt(a)-b,1,ratio);const sh=shock(mach(shockArea,'sup'));p02=p0*sh.p0Ratio;Astar2=At*shockArea/area(sh.M2);Me=mach(At*ratio/Astar2);regime='Choked · internal normal shock';
  }else{Me=mu;regime=Math.abs(b-design)<1e-6?'Choked · ideally expanded':b<design?'Choked · underexpanded exit':'Choked · overexpanded exit';}
  const sonic=state(1,p0,T0),mdot=sonic.rho*sonic.u*Astar,exit=state(Me,p02,T0),shockX=shockArea===null?null:.4+.6*(shockArea-1)/Math.max(ratio-1,1e-10);
  const A=x=>At*(x<=.4?1+1.4*(1-x/.4):1+(ratio-1)*(x-.4)/.6);
  const atX=(x,after=true)=>{const a=A(x);let st;if(!choked)st=state(mach(a/Astar),p0,T0);else if(x<=.4)st=state(mach(a/At),p0,T0);else if(shockX!==null&&(x>shockX||after&&Math.abs(x-shockX)<1e-10))st=state(mach(a/Astar2),p02,T0);else st=state(mach(a/At,'sup'),p0,T0);return {...st,x,area:a};};
  let samples=Array.from({length:121},(_,i)=>atX(i/120));
  if(shockX!==null){samples=samples.filter(v=>Math.abs(v.x-shockX)>1e-9);samples.push(atX(shockX,false),atX(shockX,true));samples.sort((a,b)=>a.x-b.x);}
  return {ok:true,regime,choked,mdot,At,ratio,exit,samples,shockX,shockArea,p02,onset,design,exitShock,backRatio:b,pressureResidual:exit.p-pb};
 }
 return {gamma,R,cp,pr,area,mach,shock,state,assess,solve};
})();
if(typeof module!=='undefined')module.exports=Compressible;
