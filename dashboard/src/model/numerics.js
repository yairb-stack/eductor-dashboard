/* Shared SI transport properties and reproducible stochastic numerics. */
const ExposureNumerics=(()=>{
 const R=287.05,VERSION='phase1-2026-11-21';
 const viscosity=T=>1.716e-5*(T/273.15)**1.5*(273.15+110.4)/(T+110.4);
 function gas(p,T){
  if(!(Number.isFinite(p)&&p>0&&Number.isFinite(T)&&T>0))throw Error('Gas properties need positive finite absolute pressure and temperature.');
  const mu=viscosity(T),rho=p/(R*T),lambda=mu/p*Math.sqrt(Math.PI*R*T/2);
  return {mu,rho,nu:mu/rho,lambda};
 }
 // WZW uses theta = 2.5 tau / Te, Te = sigma^2/epsilon (Zhou 2001, p.102).
 const lagrangianTime=(sigma,epsilon)=>epsilon>0?.4*sigma*sigma/epsilon:1;
 function normal(seed,axis,node){
  const hash=x=>{x=Math.imul(x^(x>>>16),0x21f0aaad);x=Math.imul(x^(x>>>15),0x735a2d97);return ((x^(x>>>15))>>>0)+.5;};
  const k=(seed^Math.imul(axis+1,0x9e3779b1)^Math.imul(node,0x85ebca6b))>>>0;
  return Math.sqrt(-2*Math.log(hash(k)/4294967296))*Math.cos(2*Math.PI*hash(k^0x68bc21eb)/4294967296);
 }
 // A lazy dyadic Brownian bridge. Both resolutions query the SAME W(t),
 // including adaptive wall steps, rather than drawing unrelated step noise.
 function brownian(seed,horizon,sign=1,depth=25){
  const caches=Array.from({length:3},()=>new Map());
  function value(t,axis=0){
   t=Math.max(0,Math.min(horizon,t));if(t===0)return 0;
   let lo=0,hi=horizon,a=0,b=Math.sqrt(horizon)*normal(seed,axis,1),node=1;
   if(t===horizon)return sign*b;
   const cache=caches[axis];
   for(let level=0;level<depth;level++){
    const mid=(lo+hi)/2;let c=cache.get(node);
    if(c===undefined){c=(a+b)/2+Math.sqrt((hi-lo)/4)*normal(seed,axis,node+1);cache.set(node,c);}
    if(t===mid)return sign*c;
    if(t<mid){hi=mid;b=c;node=node*2;}else{lo=mid;a=c;node=node*2+1;}
   }
   return sign*(a+(b-a)*(t-lo)/(hi-lo));
  }
  const endpoints=Array(3);
  return {value,increment:(a,b,axis=0)=>{const last=endpoints[axis],va=last&&last.t===a?last.v:value(a,axis),vb=value(b,axis);endpoints[axis]={t:b,v:vb};return vb-va;}};
 }
 function timestep({maximum,Gref,G,b,r,db,speed,remaining,TL,noisy,clearance,adaptive=false}){
  const cap=maximum*(adaptive?Math.min(20,Math.max(1,Gref/Math.max(20*G,1))):1);
  return Math.min(cap,.1*Math.min(b,db/2,r)/speed,remaining,noisy?.2*TL:Infinity,clearance>0?.2*clearance/speed:Infinity);
 }
 function firstContact(clearance,a,b,tolerance=1e-10){
  const length=Math.hypot(...b.map((v,i)=>v-a[i]));if(!length)return null;
  const at=t=>clearance(...a.map((v,i)=>v+t*(b[i]-v)));
  // Conservative sphere tracing for a signed-distance field, followed by a
  // bracketed root. It cannot jump across a wall thinner than a full step.
  let t=0,previous=0;
  for(let k=0;k<128&&t<=1;k++){
   const d=at(t);if(d<=tolerance){let lo=previous,hi=t;for(let j=0;j<35;j++){const m=(lo+hi)/2;if(at(m)>tolerance)lo=m;else hi=m;}return (lo+hi)/2;}
   previous=t;t+=Math.max(tolerance/length,.8*d/length);if(t>1&&at(1)>tolerance)return null;t=Math.min(t,1);
  }
  return null;
 }
 function wilson(fraction,n){if(!(n>0))return [0,1];const z=1.95996398454,d=1+z*z/n,c=(fraction+z*z/(2*n))/d,h=z*Math.sqrt(fraction*(1-fraction)/n+z*z/(4*n*n))/d;return [Math.max(0,c-h),Math.min(1,c+h)];}
 function precision(runs,target=.02){
  const total=runs.reduce((s,r)=>s+r.weight,0),s2=runs.reduce((s,r)=>s+(r.weight/total)**2,0),neff=1/s2;
  const wall=runs.reduce((s,r)=>s+((r.contact??r.fate==='wall')?r.weight/total:0),0),ci=wilson(wall,neff),halfWidth=(ci[1]-ci[0])/2;
  return {effectiveN:neff,halfWidth,target,met:halfWidth<=target,ci};
 }
 function finiteScalars(value,path='result',issues=[]){
  if(typeof value==='number'&&!Number.isFinite(value))issues.push(path);
  else if(value&&typeof value==='object')for(const [k,v] of Object.entries(value))if(typeof v!=='function')finiteScalars(v,path+'.'+k,issues);
  return issues;
 }
 function residence(samples){
  const exits=samples.filter(s=>s.fate==='outlet');return {definition:'Release to outlet; conditioned on reaching outlet. Wall and censored histories excluded.',count:exits.length,total:samples.length,exitTimes:exits.map(s=>s.t),exitWeights:exits.map(s=>s.powderFluxWeight??1),censored:samples.filter(s=>s.fate==='censored').length};
 }
 function wasserstein(a,b,wa=null,wb=null){
  if(!a.length||!b.length)return null;
  const pack=(v,w)=>{const total=w?w.reduce((s,x)=>s+x,0):v.length;return v.map((t,i)=>({t,w:(w?w[i]:1)/total})).sort((a,b)=>a.t-b.t);},x=pack(a,wa),y=pack(b,wb);let i=0,j=0,cx=0,cy=0,prev=Math.min(x[0].t,y[0].t),distance=0;
  while(i<x.length||j<y.length){const t=Math.min(x[i]?.t??Infinity,y[j]?.t??Infinity);distance+=Math.abs(cx-cy)*(t-prev);while(i<x.length&&x[i].t===t)cx+=x[i++].w;while(j<y.length&&y[j].t===t)cy+=y[j++].w;prev=t;}return distance;
 }

 return {VERSION,gas,viscosity,lagrangianTime,normal,brownian,timestep,firstContact,wilson,precision,finiteScalars,residence,wasserstein};
})();
if(typeof module!=='undefined')module.exports=ExposureNumerics;

'use strict';
// Gaussian isotropic well-mixed closure. The vector cross-gradient term is
// retained; using only three independent 1-D drifts is not well mixed in 3-D.
