const ExposureDispersion=(()=>{
 const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
 function gradients(field,x,y,z,a){
  const h=field.geo.dn*.035,q=[x,y,z],variance=[],sigma=[],logDensity=[],mean=[[],[],[]];
  for(let i=0;i<3;i++){const lo=q.slice(),hi=q.slice();lo[i]-=h;hi[i]+=h;const l=field.sample(...lo),r=field.sample(...hi);variance.push((r.sigma**2-l.sigma**2)/(2*h));sigma.push((r.sigma-l.sigma)/(2*h));logDensity.push((Math.log(r.rho)-Math.log(l.rho))/(2*h));['ux','uy','uz'].forEach((key,j)=>mean[j][i]=((r[key]||0)-(l[key]||0))/(2*h));}
  return {variance,sigma,logDensity,mean};
 }
 function drift(w,a,grad,slip=w){
  const s=a.sigma*a.sigma;if(!(s>1e-16))return [0,0,0];
  const conv=dot([a.ux+slip[0],a.uy+slip[1],(a.uz||0)+slip[2]],grad.variance);
  return w.map((v,i)=>.5*grad.variance[i]+.5*v*conv/s+s*grad.logDensity[i]+(grad.mean?dot(grad.mean[i],w.map((q,j)=>q-slip[j])):0));
 }
 // Limit both variance-gradient and deterministic-drift increments. This
 // is a stability guard, not a claim of second-order stochastic accuracy.
 function stepLimit(w,a,slip,grad){
  if(!(a.sigma>0))return Infinity;
  const gradient=Math.hypot(...grad.variance),d=Math.hypot(...drift(w,a,grad,slip)),ts=times(a,slip);
  return Math.min(gradient>0?.05*a.sigma/gradient:Infinity,d>0?.1*a.sigma/d:Infinity,.1*ts.parallel,.1*ts.perpendicular);
 }
 function times(a,slip){const v=Math.hypot(...slip),xi=a.sigma>0?v/a.sigma:0;return {parallel:a.TL/Math.sqrt(1+xi*xi),perpendicular:a.TL/Math.sqrt(1+4*xi*xi),direction:v?slip.map(q=>q/v):[1,0,0]};}
 function advance(w,a,slip,dt,gaussian,grad,{wellMixed=true,crossing=true}={}){
  if(!(a.sigma>0))return [0,0,0];
  const ts=crossing?times(a,slip):{parallel:a.TL,perpendicular:a.TL,direction:[1,0,0]},e=ts.direction;
  const z=[gaussian(0),gaussian(1),gaussian(2)],wp=dot(w,e),zp=dot(z,e),d=wellMixed?drift(w,a,grad,slip):[0,0,0];
  const ep=Math.exp(-dt/Math.max(ts.parallel,1e-15)),et=Math.exp(-dt/Math.max(ts.perpendicular,1e-15)),np=a.sigma*Math.sqrt(-Math.expm1(-2*dt/Math.max(ts.parallel,1e-15))),nt=a.sigma*Math.sqrt(-Math.expm1(-2*dt/Math.max(ts.perpendicular,1e-15)));
  const dp=dot(d,e),ip=-Math.expm1(-dt/Math.max(ts.parallel,1e-15))*ts.parallel,it=-Math.expm1(-dt/Math.max(ts.perpendicular,1e-15))*ts.perpendicular;
  return w.map((v,i)=>et*(v-wp*e[i])+ep*wp*e[i]+nt*(z[i]-zp*e[i])+np*zp*e[i]+(d[i]-dp*e[i])*it+dp*e[i]*ip);
 }
 // Exact change of variables for the steady isotropic Gaussian reference
 // closure: z = w/sigma. The material variance derivative cancels algebraically.
 // Carry z through zero-variance regions; never divide a finite-difference
 // variance gradient by a vanishing local variance. This is not an anisotropic
 // Reynolds-stress model. Instantaneous fluid/particle slip sets crossing time.
 function advanceNormalized(z,a,particleSlip,dt,gaussian,grad,{crossing=true,isotropic=false}={}){
  const slip=z.map((q,i)=>a.sigma*q-particleSlip[i]);
  const ts=crossing?times(a,slip):{parallel:a.TL,perpendicular:a.TL,direction:[1,0,0]};
  if(isotropic)ts.perpendicular=ts.parallel;
  const e=ts.direction,zp=dot(z,e),noise=[gaussian(0),gaussian(1),gaussian(2)],np=dot(noise,e);
  const d=z.map((_,i)=>(grad.sigma?.[i]??(a.sigma>0?grad.variance[i]/(2*a.sigma):0))+a.sigma*grad.logDensity[i]),dp=dot(d,e);
  const tp=Math.max(ts.parallel,1e-15),tt=Math.max(ts.perpendicular,1e-15),ep=Math.exp(-dt/tp),et=Math.exp(-dt/tt),ip=-Math.expm1(-dt/tp)*tp,it=-Math.expm1(-dt/tt)*tt;
  return z.map((v,i)=>et*(v-zp*e[i])+ep*zp*e[i]+Math.sqrt(-Math.expm1(-2*dt/tt))*(noise[i]-np*e[i])+Math.sqrt(-Math.expm1(-2*dt/tp))*np*e[i]+it*(d[i]-dp*e[i])+ip*dp*e[i]);
 }
 return {gradients,drift,times,stepLimit,advance,advanceNormalized};
})();
if(typeof module!=='undefined')module.exports=ExposureDispersion;

'use strict';
// Fixed-pivot sectional aggregation (Kumar & Ramkrishna, 1996).
// q_i is number per kg of conserved powder, C is local powder kg / m³.
// Thus dilution changes rates, never the parcel's mass. Unbounded growth is
// recorded in an overflow reservoir rather than silently discarded.
