/* Compressible reduced eductor estimate, SI units.
 * Isentropic inlet restrictions + constant-area mixing (mass, axial momentum,
 * total energy), with subsonic area change or critical mixing and an internal normal shock.
 * Powder is pressureless, co-moving inertia with no thermal exchange here.
 * Geometry-specific mixing/separation efficiencies remain uncalibrated.
 */
const EductorNetwork=(()=>{
 const R=287.05,gamma=1.4,cp=gamma*R/(gamma-1),crit=(2/(gamma+1))**(gamma/(gamma-1));
 const area=d=>Math.PI*d*d/4,clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
 function inlet(p0,p,T0,A,md=0){
  if(!(p0>p&&p>0&&T0>0&&A>0))return {m:0,u:0,T:T0,p:p0,p0,T0,A,choked:false,M:0,md};
  const pe=Math.max(p,p0*crit),T=T0*(pe/p0)**((gamma-1)/gamma),ug=Math.sqrt(2*cp*(T0-T)),m0=A*pe/(R*T)*ug;
  const m=md>0?2*m0*m0/(Math.hypot(md,2*m0)+md):m0,z=m>0?md/m:0,u=ug/Math.sqrt(1+z);
  return {m,u,T,p:pe,p0,T0,A,md,choked:p/p0<=crit,M:u/Math.sqrt(gamma*R*T/(1+z))};
 }
 function flux(M,gamma=1.4){return M*(1+(gamma-1)*M*M/2)**(-(gamma+1)/(2*(gamma-1)));}
 function duct(m,md,T0,p0,A,branch='subsonic'){
  const z=md/m,scale=A*p0*Math.sqrt(gamma/(R*T0))/Math.sqrt(1+z),target=m/scale;
  if(!(target>=0)||target>flux(1,gamma)*(1+1e-9))return null;
  const supersonic=branch==='supersonic';let lo=supersonic?1:0,hi=supersonic?50:1;if(supersonic&&target<flux(hi,gamma))return null;for(let i=0;i<(supersonic?70:60);i++){const mid=(lo+hi)/2;if(supersonic?flux(mid,gamma)>target:flux(mid,gamma)<target)lo=mid;else hi=mid;}
  const M=(lo+hi)/2,T=T0/(1+(gamma-1)*M*M/2),p=p0*(T/T0)**(gamma/(gamma-1)),u=M*Math.sqrt(gamma*R*T/(1+z));
  return {m,md,T0,p0,A,M,T,p,u,rho:p/(R*T)};
 }
 // Ideal normal shock for the existing homogeneous, co-moving powder closure.
 // Constant solids loading scales R and cp equally, so gamma is unchanged.
 function normalShock(up){
  if(!up||up.M<1-1e-8)return null;
  const M2=up.M*up.M,ratio=(gamma+1)*M2/(2+(gamma-1)*M2),p=up.p*(2*gamma*M2-(gamma-1))/(gamma+1),T=up.T*(p/up.p)/ratio,u=up.u/ratio,
   M=u/Math.sqrt(gamma*R*T/(1+up.md/up.m)),p0=p*(up.T0/T)**(gamma/(gamma-1)),rho=p/(R*T),mt=up.m+up.md;
  return {...up,p,T,u,M,p0,rho,residuals:{mass:up.m-rho*up.A*u,momentum:p*up.A+mt*u-up.p*up.A-mt*up.u,energy:up.m*cp*(T-up.T)+.5*mt*(u*u-up.u*up.u)}};
 }
 function shockDiffuser(mix,Ae,backPressure){
  if(!(Ae>mix.A))return null;
  const at=A=>{const before=duct(mix.m,mix.md,mix.T0,mix.p0,A,'supersonic'),after=normalShock(before);if(!after)return null;const out=duct(mix.m,mix.md,mix.T0,after.p0,Ae);return out?{before,after,out,area:A,error:out.p-backPressure}:null;};
  let lo=mix.A,hi=Ae,a=at(lo),b=at(hi);if(!a||!b||a.error < -1e-6||b.error > 1e-6)return null;
  let result;
  for(let i=0;i<75;i++){const mid=(lo+hi)/2;result=at(mid);if(!result)return null;if(Math.abs(result.error)<1e-7)break;if(result.error>0)lo=mid;else hi=mid;}
  return result&&Math.abs(result.error)<.02?result:null;
 }
 function mixing(P,A,primary,secondary,bleed,md,angle,K,diagnostics=null){
  const fail=constraint=>{if(diagnostics)diagnostics.constraint=constraint;return null;};
  const m=primary.m+secondary.m+bleed.m,mt=m+md;
  if(!(m>0&&A>primary.A))return fail('mixing-area-or-mass');
  const T0=(primary.m*primary.T0+secondary.m*secondary.T0+bleed.m*bleed.T0)/m;
  const I=P*A+primary.m*primary.u+(primary.p-P)*primary.A+((secondary.m+md)*secondary.u+(secondary.p-P)*secondary.A)*Math.cos(angle*Math.PI/180);
  const B=mt*(1-R/(2*cp)+K/2),C=m*R*T0,disc=I*I-4*B*C;
  if(!(disc>=0&&I>0))return fail('momentum-discriminant');
  const u=2*C/(I+Math.sqrt(disc)),T=T0-mt/m*u*u/(2*cp),p=m*R*T/(A*u),rho=p/(R*T),M=u/Math.sqrt(gamma*R*T/(1+md/m));
  if(!(p>0&&T>0&&M<1+1e-8))return fail('mixing-choking-or-temperature');
  const p0=p*(T0/T)**(gamma/(gamma-1)),wall=K/2*mt*u;
  const entropy=m*(cp*Math.log(T0)-R*Math.log(p0))-[primary,secondary,bleed].reduce((s,v)=>s+v.m*(cp*Math.log(v.T0)-R*Math.log(v.p0)),0);
  if(entropy < -1e-7)return fail('entropy');
  return {m,md,T0,p0,A,u,T,p,rho,M,wall,entropy,
   residuals:{mass:m-rho*A*u,momentum:p*A+mt*u+wall-I,energy:m*cp*T+.5*mt*u*u-m*cp*T0}};
 }
 function dimensions(g,loss=1){
  const ap=(g.nports||1)*area(g.db/1000),af=area((g.df||g.db)/1000),british=g.family==='british';
  let As=british?ap/Math.sqrt(1+1.3*loss):1/Math.sqrt(2*loss/(af*af)+(1+.3*loss)/(ap*ap));
  if(g.key==='revA')As=ap/Math.sqrt(1+.31*loss);
  if(g.key==='revB1')As=ap/Math.sqrt(1+.20*loss);
  if(g.key==='revB2')As=ap/Math.sqrt(1+.28*loss);
  // The insert constrains the secondary passage; the primary jet has its own nozzle.
  if(british&&g.insert){const annulus=Math.PI*((g.rt/1000)**2-(g.dn/2000)**2);if(!(annulus>0))return null;As=1/Math.sqrt(1/(As*As)+(1+.59*loss)/(annulus*annulus));}
  return {A:area(g.d/1000),Ae:area(g.de/1000),As,Ah:british?(g.nh||8)*area((g.dh||8)/1000)/Math.sqrt(1+1.5*loss):0,
   An:area(g.dn/1000)*(british?.94:.82),K:.02*loss*g.L/g.d+(g.de<g.d?.2*loss:0)};
 }
 function stage(g,P,p0,Ps,Ts,Po,Tm,md,loss,secondary=null){
  const d=dimensions(g,loss);if(!d)return null;
  const n=p0===null?{m:0,u:0,T:Tm,T0:Tm,p:Ps,p0:Ps,A:0,M:0,choked:false}:inlet(p0,P,Tm,d.An),s=secondary||inlet(Ps,P,Ts,d.As,md),h=d.Ah?inlet(Po,P,293.15,d.Ah):{m:0,u:0,T0:293.15,p0:Po,T:293.15,M:0,choked:false};
  if(!((n.m>0||p0===null)&&s.m>1e-10))return null;
  const mix=mixing(P,d.A,n,s,h,md,g.angle??(g.family==='british'?90:45),d.K);if(!mix)return null;
  const out=duct(mix.m,md,mix.T0,mix.p0,d.Ae);if(!out)return null;
  return {g,P,n,s,h,mix,out,d};
 }
 // Below the fully subsonic branch's backpressure range, let the mixed
 // stream become sonic at the nominal bore and solve a normal shock in its
 // expanding diffuser. This is a separate critical operating branch, not an
 // imposed carrier flow or a change to the user's boundary pressures.
 function criticalStage(g,p0,Ps,Ts,Po,Tm,md,loss){
  const d=dimensions(g,loss);if(!d||!(d.Ae>d.A))return null;
  const critical=root(P=>{
   const n=inlet(p0,P,Tm,d.An),s=inlet(Ps,P,Ts,d.As,md),h=d.Ah?inlet(Po,P,293.15,d.Ah):{m:0,u:0,T0:293.15,p0:Po,T:293.15,M:0,choked:false};
   if(!n.choked||!(s.m>1e-10))return null;
   const m=n.m+s.m+h.m,mt=m+md,T0=(n.m*n.T0+s.m*s.T0+h.m*h.T0)/m,T=T0/(1+(gamma-1)/2),u=Math.sqrt(gamma*R*T/(1+md/m)),pressure=m*R*T/(d.A*u),totalP=pressure*(T0/T)**(gamma/(gamma-1)),rho=pressure/(R*T),wall=d.K/2*mt*u;
   const I=P*d.A+n.m*n.u+(n.p-P)*n.A+((s.m+md)*s.u+(s.p-P)*s.A)*Math.cos((g.angle??90)*Math.PI/180),residual=pressure*d.A+mt*u+wall-I;
   const entropy=m*(cp*Math.log(T0)-R*Math.log(totalP))-[n,s,h].reduce((v,q)=>v+q.m*(cp*Math.log(q.T0)-R*Math.log(q.p0)),0);
   if(!(pressure>0&&T>0)||entropy < -1e-7)return null;
   const mix={m,md,T0,p0:totalP,A:d.A,u,T,p:pressure,rho,M:1,wall,entropy,residuals:{mass:m-rho*d.A*u,momentum:residual,energy:m*cp*T+.5*mt*u*u-m*cp*T0}};
   return {g,P,n,s,h,mix,d,error:residual/d.A};
  },Math.min(Ps,p0),{dense:true,accept:q=>!!shockDiffuser(q.mix,d.Ae,Po)});
  if(!critical)return null;
  const shock=shockDiffuser(critical.mix,d.Ae,Po);if(!shock)return null;
  return {...critical,out:shock.out,shock,error:shock.error};
 }
 function root(evaluate,limit,{seeds=[],dense=false,accept=()=>true}={}){
  // A root can sit between a choking/entropy boundary and the first coarse
  // admissible sample. Refine BOTH edges of every detected admissible interval
  // before declaring it unbracketed. Never interpolate through an invalid state.
  const ratios=[...seeds.flatMap(v=>[v/limit,...[.00001,.0001,.001,.01,.05].flatMap(d=>[v/limit-d,v/limit+d])]).filter(v=>v>0&&v<1),...(dense?Array.from({length:1024},(_,i)=>(i+.5)/1024):[]),...Array.from({length:90},(_,i)=>.001+(i/89)*.979),...Array.from({length:65},(_,i)=>1-10**(-1.7-i*.1))].sort((a,b)=>a-b);
  let prev=null,best=null;
  const tolerance=Math.max(.02,limit*1e-7);
  const sample=x=>{const value=evaluate(x),r=value&&Number.isFinite(value.error)?value:null;if(r&&(!best||Math.abs(r.error)<Math.abs(best.error)))best=r;return {x,r};};
  const boundary=(good,bad)=>{for(let k=0;k<42;k++){const mid=sample((good.x+bad.x)/2);if(mid.r)good=mid;else bad=mid;}return good;};
  const bracket=(a,b)=>{
   if(!a.r||!b.r)return null;
   if(Math.abs(a.r.error)<1e-7)return a.r;if(Math.abs(b.r.error)<1e-7)return b.r;
   if(Math.sign(a.r.error)===Math.sign(b.r.error))return null;
   let answer=null;
   for(let k=0;k<65;k++){const mid=sample((a.x+b.x)/2);if(!mid.r)return null;answer=mid.r;if(Math.abs(answer.error)<1e-7)break;if(Math.sign(answer.error)===Math.sign(a.r.error))a=mid;else b=mid;}
   return answer&&Math.abs(answer.error)<tolerance?answer:null;
  };
  for(const v of ratios){const point=sample(limit*v);
   if(point.r&&Math.abs(point.r.error)<1e-7&&accept(point.r))return point.r;
   if(prev){let found=null;
    if(prev.r&&point.r)found=bracket(prev,point);
    else if(!prev.r&&point.r)found=bracket(boundary(point,prev),point);
    else if(prev.r&&!point.r)found=bracket(prev,boundary(prev,point));
    if(found&&accept(found))return found;
   }
   prev=point;
  }
  return best&&Math.abs(best.error)<tolerance&&accept(best)?best:null;
 }
 function secondaryPressure(m,md,p0,T0,A){
  const cap=inlet(p0,p0*crit,T0,A,md);if(m>cap.m*(1+1e-9))return null;
  let lo=p0*crit,hi=p0;for(let i=0;i<60;i++){const p=(lo+hi)/2;if(inlet(p0,p,T0,A,md).m>m)lo=p;else hi=p;}return (lo+hi)/2;
 }
 const trainGeometries=[{key:'revB1',dn:3.5,db:34.8,d:34.8,de:34.8,L:110.3,angle:90,nports:1,family:'branched'},
  {key:'revB2',dn:3.5,db:22.1,d:22.1,de:22.1,L:61.7,angle:45,nports:1,family:'branched'}];
 function solve(g,p){
  const po=p.ps*1e5,ps=po,Tm=p.temp+273.15,p0=p.p0*1e5,md=p.feed/3600,loss=p.loss??1;
  if(![po,p0,Tm,md,loss].every(Number.isFinite)||!(po>0&&p0>po&&Tm>0&&md>=0&&loss>0))return {ok:false,reason:'Compressible estimate needs positive inputs and motive pressure above the shared vessel / outlet pressure.'};
  let result;
  if(g.key?.startsWith('revB')){
   const gs=trainGeometries.map(v=>({...v}));gs[1].L=p.ext!=null?61.7+p.ext*1000:gs[1].L;
   result=root(P=>{const a=stage(gs[0],P,p.off1?null:(p.p01??p.p0)*1e5,ps-(p.feedDrop||0),293.15,po,Tm,md,loss);if(!a)return null;
    const pt=a.out.p0-(p.interDrop||0),P2=secondaryPressure(a.mix.m,md,pt,a.mix.T0,dimensions(gs[1],loss).As);if(!P2)return null;
    const incoming=inlet(pt,P2,a.mix.T0,dimensions(gs[1],loss).As,md),b=stage(gs[1],P2,p0,pt,a.mix.T0,po,Tm,md,loss,incoming);if(!b)return null;
    return {stages:[a,b],error:b.out.p-po};},Math.min(ps-(p.feedDrop||0),(p.p01??p.p0)*1e5));
  }else result=root(P=>{const a=stage(g,P,p0,ps-(p.feedDrop||0),293.15,po,Tm,md,loss);return a?{stages:[a],error:a.out.p-po}:null;},Math.min(ps-(p.feedDrop||0),p0));
  if(!result&&!g.key?.startsWith('revB')){const a=criticalStage(g,p0,ps-(p.feedDrop||0),293.15,po,Tm,md,loss);if(a)result={stages:[a],error:a.error};}
  if(!result)return {ok:false,reason:'The compressible model found no outlet-matching branch after refining the admissible pressure boundaries. No artificial airflow has been substituted.',attempted:true};
  const stages=result.stages.map(a=>toFlow(a,p,po,result.error));
  return {...stages[g.key==='revB2'?1:0],stages};
 }
 function toFlow(a,p,po,residual){
  const {g,P,n,s,h,mix,out}=a,md=p.feed/3600;
  const power=n.m*cp*n.T0*(1-(P/n.p0)**((gamma-1)/gamma))+s.m*cp*s.T0*(1-(P/s.p0)**((gamma-1)/gamma))+h.m*cp*h.T0*(1-(P/h.p0)**((gamma-1)/gamma));
  return {ok:true,manual:false,compressible:true,closure:a.shock?'compressible-choked-diffuser':'compressible-mixing',caution:true,reason:a.shock?'Automatic choked-mixing / internal-shock estimate; uncalibrated homogeneous flow.':'Automatic compressible mixing estimate; uncalibrated complete-mixing closure.',
   mp:n.m,ms:s.m,mh:h.m,mg:mix.m,md,p0:n.p0,ps:po,po:out.p,systemOutletP:po,inletP:s.p0,Tin:s.T0,T0:n.T0,T:out.T,portP:P,throatP:mix.p,
   totalT:mix.T0,uJet:n.m>0?Math.sqrt(2*cp*n.T0*(1-(P/n.p0)**((gamma-1)/gamma))):0,power,choked:n.choked,
   network:{regime:a.shock?'Choked mixing · internal diffuser shock':'Subsonic mixing / outlet',mixingMach:mix.M,outletMach:out.M,shock:a.shock?{area:a.shock.area,diameter:Math.sqrt(4*a.shock.area/Math.PI),before:a.shock.before,after:a.shock.after}:null,bulkMach:Math.max(mix.M,out.M,s.M,h.M,a.shock?.before.M||0),motiveThrust:n.m*n.u+(n.p-P)*n.A,motiveChoked:n.choked,secondaryChoked:s.choked,bleedChoked:h.choked,
    pressureResidual:residual,residuals:mix.residuals,entropyProduction:mix.entropy,mixingPressure:mix.p,mixingTemperature:mix.T,mixingVelocity:mix.u,outletVelocity:out.u,secondaryVelocity:s.u,bleedVelocity:h.u,
    assumptions:a.shock?'Critical complete mixing at the nominal bore (M = 1), with lumped wall loss upstream; ideal expansion and one normal shock in the diffuser matched to the outlet pressure. Powder is homogeneous co-moving inertia, without cp exchange. Shock area is estimated; shock cells and separation are unresolved. Particle histories still use the smooth analytical pressure field, without resolved shock-particle dynamics.':'Complete mixing at the nominal bore; lumped bleed entry; ideal subsonic outlet area change. Powder adds co-moving inertia, without thermal exchange. Insert secondary area is restricted, but shock cells, jet blockage and separation are not resolved.'}};
 }
 function assess(r,g,p){
  if(!r?.supported)return {needed:true,reason:'Low-Mach screening did not pass: '+(r?.reason||'reduced balance unavailable')};
  if(r.bulkMach>=.3)return {needed:true,reason:'Estimated bulk Mach ≥ 0.30'};
  if(g&&p){
   const portArea=(g.nports||1)*area(g.db/1000),feedArea=area((g.df||g.db)/1000);
   const restricted=g.insert?Math.PI*((g.rt/1000)**2-(g.dn/2000)**2):Infinity;
   const A=Math.min(portArea,feedArea,restricted),rho=p.ps*1e5/(R*293.15),M=r.ms/(rho*A)*Math.sqrt(1+(r.md||0)/r.ms)/Math.sqrt(gamma*R*293.15);
   if(M>=.3)return {needed:true,reason:'Estimated feeder / secondary restriction Mach ≥ 0.30'};
  }
  if(r.caution)return {needed:true,reason:'Pressure / density variation exceeds the reduced-flow screen'};
  return {needed:false,reason:'Low-Mach balance passes its screening checks'};
 }
 function legacyResult(f,g,base,po){
  if(!f.ok)return {...base,supported:false,compressibleAttempted:true,reason:f.reason};
  const d=g.d/1000,L=g.L/1000,st={...base.st,d,db:g.db/1000*Math.sqrt(g.nports||1),L,Ts:f.Tin};
  const jetT=f.T0-f.uJet*f.uJet/(2*cp);
  const n={...base.n,mdot:f.mp,p0:f.p0,Tm:f.T0,choked:f.choked,uj:f.uJet,mj:f.uJet/Math.sqrt(gamma*R*jetT),power:f.power,uEq:f.mp>0?f.network.motiveThrust/f.mp:0,fNet:f.network.motiveThrust,dEq:f.mp>0?Math.sqrt(4*f.mp/(Math.PI*f.portP/(R*jetT)*f.uJet)):0,aEff:f.mp>0?area(g.dn/1000)*(g.family==='british'?.94:.82):0};
  n.S=n.fNet+f.portP*n.aEff;
  const r={...base,st,n,mp:f.mp,ms:f.ms,mh:f.mh,mg:f.mg,md:f.md,Tmix:f.T,omega:f.ms/f.mp,p1:f.portP-po,pCh:f.portP-po,pGap:f.portP-po,pThroat:f.throatP-po,pMix:f.throatP-po,pe:f.po-po,pt:f.po*(f.totalT/f.T)**(gamma/(gamma-1))-po,
   vb:f.network.secondaryVelocity,vPort:f.network.secondaryVelocity,vHole:f.network.bleedVelocity,vExit:f.network.outletVelocity,vMix:f.network.mixingVelocity,vThroat:f.network.mixingVelocity,vJet:f.network.outletVelocity,vAnn1:f.network.secondaryVelocity,
   supported:true,caution:f.modelVersion===2?f.compressible:true,reason:f.reason,compressible:f.compressible,modelVersion:f.modelVersion,bulkMach:f.network.bulkMach,network:f.network,flow:f,beta:1,loading:f.md/f.mg,loadFeed:f.md/f.ms,feedLoading:f.md/f.ms,conc:f.md/f.mg*f.po/(R*f.T),ks:base.ks||0};
  r.lAttach=d/2/.12;r.lEnt=f.mg/(.12*Math.sqrt(Math.PI*f.po/(R*f.T)*n.fNet));r.eps=Math.max(0,f.power-.5*(f.mg+f.md)*r.vExit**2)/(f.po/(R*f.T)*area(d)*L);r.shear=Math.sqrt(r.eps/(1.81e-5/(f.po/(R*f.T))));r.tRes=L/r.vExit;
  if(base.rows){r.portAir=f.ms/(g.nports||1);r.portPowder=f.md/(g.nports||1);r.rows=base.rows.map((row,i)=>{let air=i===0?f.mp:i===1?f.ms:i===2?r.portAir:i===3?f.mp+f.ms:f.mg,powder=i===0?0:i===2?r.portPowder:f.md,pressure=i===0?f.p0:i===1?po:i<4?f.portP:i===base.rows.length-1?f.po:f.throatP;return {...row,air,powder,pressure,ratio:powder/air,fraction:powder/(air+powder)};});}
  if(f.circuit?.flowUsed&&f.circuit.branches.length>1&&r.rows){r.rows[2].name='Mean suction branch';r.rows[2].detail='Arithmetic mean only; individual CAD branch flows are listed separately.';}
  return r;
 }
 return {inlet,duct,normalShock,shockDiffuser,mixing,criticalStage,dimensions,solve,solvePressureRoot:root,assess,legacyResult,R,gamma,cp};
})();
if(typeof module!=='undefined')module.exports=EductorNetwork;

const EDUCTOR_BUILD={"codeSHA256": "5511e2bda37281b3b16892a4b15860aa5a00bccf1a791e28d576fdcae09966c5", "sources": {"exposure/assemblies.js": "74ddd6abfb9dafeb274ec7a229dcb01f127042d7978735d85e8d17270a2a0cb0", "exposure/collisions.js": "480783928f9f4938e64eeab427f232b1734103ea573e396182276c736d7597a6", "exposure/dispersion.js": "401eb2de5eb540d48c2cfb9817801b29f5e6b8c0b0c35997914085768c737bcb", "exposure/geometry.json": "bc79d875e51bdf344b55f230a278079f95da278221422e7d094cc55dbb3cb009", "exposure/model.js": "c480474d489eeaf9756685742b2f8c9868caa209dfa8ad901798d80a1d992abd", "exposure/network.js": "fdd316d6679443c2975eddcdf7886bea05f56c03438c6faee76c99358af9e162", "exposure/numerics.js": "c366f11a66e0ebd24e2b1ed822ca83e769dfcb4959c79d2bccb958014575f8ca", "exposure/population.js": "c29e6414fa916fbaaf1e363f5a8fdcaccf27f5e5e95b38f1db7b5044e7ddce2c", "exposure/runtime.js": "0407d9bcbd438591bfc83dfc980b473a43b9f03bc948c2b54abaa76427bbbb3c", "exposure/worker.js": "646cd58ce19f90a74b79e8e92ead1dc22ee069f3718c325503a16f23895f25e6", "phases/calibration-candidate.json": "2d55517d65f11761b6e685c24170a72a4b661edb0faf2e711889e968251e8252", "phases/walls.json": "87ab6bd67451bedeb826bef9e6f46f7a5ed92f2b7cc90521c18b363f4c8b9a2b", "upgrades/calibration.json": "062a6020a1844932c81a0f115c5898ebfb1e37b8e76e8d36ee1ff5fdc8bb7d79", "upgrades/continuation.js": "08dd4234d74e7adc89595d514d7d2a359dea1b8abd82719168d3c1d96d375d7b", "upgrades/jet.js": "653429c765d91545f140859069bbc8a8dbc0aff6baa48a28311c40f5207e4dd7", "upgrades/physics.js": "4426523bdb8a5b11a502457c84d50d530d9901decb743be64c87dc4ce74ebf05", "upgrades/spatial.js": "e041ed4a4a3c667cb1e472cb37265321420566626da41234392c7dc0bdb54f25", "upgrades/topology.json": "ac17cdfc828c33a27681e63ecb03a0fc8b16bdc101568103896d40baeb6f80b7", "upgrades/walls.js": "16bc58c14373ac9e18e9da658af6c2567d5f5efd3f5664db79f2c7f381d81f11"}, "buildDate": "2026-11-23T03:00:45.024772+00:00"};const EDUCTOR_WALLS={};const EDUCTOR_CANDIDATE={
  "method": "Family-level log-space MAP calibration with additive log discrepancy and Laplace prediction intervals; NOT a sampled Bayesian posterior.",
  "priors": {
    "logFeederSD": 1.5,
    "logMixingSD": 1.5,
    "logDiscrepancySD": 0.1,
    "observationLogSD": 0.05
  },
  "folds": [
    {
      "id": "british35",
      "family": "british",
      "sourceHash": "cfa61df0b09ba35b0b9be170abb80332416b2869def9b1cf6375ce97fea45a83",
      "training": [
        "british25"
      ],
      "parameters": [
        0.8819595252893219,
        0.027996469716078387,
        -0.027734375000000012,
        0.030761718750000007
      ],
      "raw": [
        0.0058742264882276395,
        0.018793895674413825
      ],
      "prediction": [
        0.00571354696047318,
        0.019381012267633273
      ],
      "targets": [
        0.0057894359934826,
        0.02011617896767593
      ],
      "errorsPct": [
        -1.3108191038790462,
        -3.6546040936699553
      ],
      "intervals": [
        [
          0.004936725271004886,
          0.006612605943714465
        ],
        [
          0.016862691474778202,
          0.02227542602436733
        ]
      ],
      "covered": [
        true,
        true
      ],
      "parameterCorrelation": -0.7317921863771872
    },
    {
      "id": "british25",
      "family": "british",
      "sourceHash": "fc554dcd6db54b484064673d0265a10308083cb349658c945df0f5dd2127c010",
      "training": [
        "british35"
      ],
      "parameters": [
        0.9798110877893225,
        -0.09016759278392161,
        -0.02646484375000001,
        0.036328125
      ],
      "raw": [
        0.005342882969322652,
        0.013607037586503057
      ],
      "prediction": [
        0.005203339054790178,
        0.014110444297480672
      ],
      "targets": [
        0.005086637371566818,
        0.013839213085114885
      ],
      "errorsPct": [
        2.2942795937390237,
        1.9598745296979025
      ],
      "intervals": [
        [
          0.004542703426654457,
          0.0059600495071377236
        ],
        [
          0.012343645440348966,
          0.016130132644726704
        ]
      ],
      "covered": [
        true,
        true
      ],
      "parameterCorrelation": -0.7621404778898444
    },
    {
      "id": "cad45-sst-coarse-v3",
      "family": "branched",
      "sourceHash": "927cd55e0c520552deacee73707555e03b992b00ebea1f8af09b843eb88e0c35",
      "training": [
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.444140625,
        0.1401376986857903,
        -0.025439453124999952,
        0.004687500000000001
      ],
      "raw": [
        0.01595128757662176,
        0.018801659691344844
      ],
      "prediction": [
        0.01555061361556474,
        0.018889999355481146
      ],
      "targets": [
        0.0169282678493279,
        0.020094052970119774
      ],
      "errorsPct": [
        -8.138187828932875,
        -5.992089383008381
      ],
      "intervals": [
        [
          0.013992567884499684,
          0.01728214476546983
        ],
        [
          0.0169916189639062,
          0.021000475376011254
        ]
      ],
      "covered": [
        true,
        true
      ],
      "parameterCorrelation": 0.8075397797885335
    },
    {
      "id": "cad45-p2bar-sst-v1",
      "family": "branched",
      "sourceHash": "9c1786665923443c5c31711d70063c7421223403d506481c75daca88ef9a916e",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.36035156249999983,
        -0.03906152006420971,
        -0.05146484374999996,
        -0.015478515625000005
      ],
      "raw": [
        0.012534779965742222,
        0.014435028042224279
      ],
      "prediction": [
        0.011905998325831064,
        0.014213315551803355
      ],
      "targets": [
        0.011817516366821166,
        0.013923453920772641
      ],
      "errorsPct": [
        0.7487356586898475,
        2.0818227480055462
      ],
      "intervals": [
        [
          0.010714918327135234,
          0.013229479852936168
        ],
        [
          0.012782367967039299,
          0.01580445340769887
        ]
      ],
      "covered": [
        true,
        true
      ],
      "parameterCorrelation": 0.8060136028137729
    },
    {
      "id": "cad45-p5bar-sst-v1",
      "family": "branched",
      "sourceHash": "7cd0fa1025054b864cbe705720e92083fc102de4d778642bdeb2adeb7cc38c1a",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.5121093750000005,
        0.2739267611857903,
        -0.003027343750000109,
        0.022460937500000052
      ],
      "raw": [
        0.02041822991462932,
        0.025168850105834462
      ],
      "prediction": [
        0.020356510384165912,
        0.025740562639505693
      ],
      "targets": [
        0.023706507914533916,
        0.02900162824184617
      ],
      "errorsPct": [
        -14.131130331153486,
        -11.244422468787862
      ],
      "intervals": [
        [
          0.01830385471738099,
          0.022639357742888996
        ],
        [
          0.023159054282688652,
          0.028609828230058212
        ]
      ],
      "covered": [
        false,
        false
      ],
      "parameterCorrelation": 0.8083168219043618
    },
    {
      "id": "cad45-ds10mm-p3bar-sst-inlet-v2",
      "family": "branched",
      "sourceHash": "0aceedafba7ee3c23b2a8fbe0e97f2623875f68c83d7af08af1b9bf2e72fee8a",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.4107421874999999,
        -0.09814355131420971,
        -0.06162109374999997,
        -0.027099609375000028
      ],
      "raw": [
        0.0053537080980751205,
        0.008204080212798207
      ],
      "prediction": [
        0.005033765586194226,
        0.007984738307088488
      ],
      "targets": [
        0.005253446712289514,
        0.008397425346610617
      ],
      "errorsPct": [
        -4.181657074419021,
        -4.914447256011611
      ],
      "intervals": [
        [
          0.004449665404008581,
          0.005694539628513703
        ],
        [
          0.007124301295611592,
          0.008949094540956406
        ]
      ],
      "covered": [
        true,
        true
      ],
      "parameterCorrelation": 0.7161137077628039
    },
    {
      "id": "cad45-ds40mm-p3bar-sst-v1",
      "family": "branched",
      "sourceHash": "c9ea0563f02814b0c96aec7fbadbb0a4f18098eb5b5b7b5676ad12234a8f530d",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.2150390625,
        -0.01171777006420971,
        -0.073486328125,
        -0.03364257812499999
      ],
      "raw": [
        0.016812196031755756,
        0.01966256814647884
      ],
      "prediction": [
        0.015621032590582031,
        0.019012072194615787
      ],
      "targets": [
        0.01784467558608473,
        0.020988287002385555
      ],
      "errorsPct": [
        -12.461100706345674,
        -9.41579847628894
      ],
      "intervals": [
        [
          0.014020500403531568,
          0.017404276036721313
        ],
        [
          0.017090845213557516,
          0.02114926936712008
        ]
      ],
      "covered": [
        false,
        true
      ],
      "parameterCorrelation": 0.78930227449065
    },
    {
      "id": "cad45-dm10mm-p3bar-sst-outlet-v2",
      "family": "branched",
      "sourceHash": "7c9e64af7b20af363863347d040fe38ab16b25d8cbc1281962c70411fe324abf",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.43945312499999983,
        -0.4316396450642097,
        -0.06455078124999986,
        -0.03813476562500002
      ],
      "raw": [
        0.006139966080001176,
        0.008990338194724262
      ],
      "prediction": [
        0.005756147628551414,
        0.008653948589670009
      ],
      "targets": [
        0.004243058571968738,
        0.007407365214469467
      ],
      "errorsPct": [
        35.66033866651568,
        16.828971423813144
      ],
      "intervals": [
        [
          0.005152789437627529,
          0.0064301551465944275
        ],
        [
          0.007777265272075454,
          0.009629455029848305
        ]
      ],
      "covered": [
        false,
        false
      ],
      "parameterCorrelation": 0.7711785176714785
    },
    {
      "id": "cad45-dm40mm-p3bar-sst-long-v2",
      "family": "branched",
      "sourceHash": "2979b57a7c29098918dc99c2ea4b7f133e6756c5167b6a6e1e627723e0c5427d",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.20039062500000004,
        -0.06249902006420971,
        -0.06235351562499999,
        -0.02426757812500003
      ],
      "raw": [
        0.01673093941800278,
        0.019581311532725864
      ],
      "prediction": [
        0.015719565544502354,
        0.019111840033150338
      ],
      "targets": [
        0.014912497208186094,
        0.01808412526563588
      ],
      "errorsPct": [
        5.412026738708975,
        5.68296642728614
      ],
      "intervals": [
        [
          0.013885679008252735,
          0.01779565412401099
        ],
        [
          0.016939536043502793,
          0.02156271745074301
        ]
      ],
      "covered": [
        true,
        true
      ],
      "parameterCorrelation": 0.7564467590788049
    },
    {
      "id": "cad90-angle-only-sst-v1",
      "family": "branched",
      "sourceHash": "8a11b5e4385073447ad1e0453f7a42754943668338d73421f51f99421c13d8a0",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.33867187500000007,
        -0.05566308256420971,
        -0.06064453125000002,
        -0.023730468750000042
      ],
      "raw": [
        0.013374668199193427,
        0.016225040313916514
      ],
      "prediction": [
        0.012587672395083486,
        0.0158445450232701
      ],
      "targets": [
        0.013087960487193136,
        0.016229029400963223
      ],
      "errorsPct": [
        -3.822506131487735,
        -2.3691150480650736
      ],
      "intervals": [
        [
          0.011331262761635442,
          0.0139833926420287
        ],
        [
          0.014273293247391305,
          0.01758876544068181
        ]
      ],
      "covered": [
        true,
        true
      ],
      "parameterCorrelation": 0.8130229020323216
    },
    {
      "id": "cad90-d38p1-dj2-p3g-open-sst-v2",
      "family": "branched",
      "sourceHash": "7efd7c8f9eaa33efecbe55cb85aa2d613df04163843414c5850ad3d1389c75df",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1"
      ],
      "parameters": [
        0.5384765625000003,
        0.2163095736857903,
        0.018847656250000018,
        0.04262695312500001
      ],
      "raw": [
        0.02226467013013308,
        0.024696987668030114
      ],
      "prediction": [
        0.022688286527205902,
        0.02577250516927548
      ],
      "targets": [
        0.019971019018703552,
        0.022623577870895088
      ],
      "errorsPct": [
        13.606053381440052,
        13.918785597707961
      ],
      "intervals": [
        [
          0.020394876249960357,
          0.025239591514638064
        ],
        [
          0.023225754450788408,
          0.028598512229504092
        ]
      ],
      "covered": [
        false,
        false
      ],
      "parameterCorrelation": 0.8180946066007959
    }
  ],
  "fits": {
    "british": {
      "parameters": [
        1.0419204627893222,
        -0.09407384278392161,
        -0.028124999999999983,
        0.02753906249999997
      ],
      "feederResistance": 2.834655641019057,
      "mixingLoss": 0.9102155469091614,
      "logDiscrepancy": [
        -0.028124999999999983,
        0.02753906249999997
      ],
      "training": [
        "british35",
        "british25"
      ],
      "covariance": [
        [
          0.763910677593919,
          -0.3247924572811161,
          0.04124491360514729,
          -0.04752564184773692
        ],
        [
          -0.32479245728111644,
          0.23910022216532806,
          0.005904300189615462,
          0.03686627308173652
        ],
        [
          0.041244913605147177,
          0.005904300189615533,
          0.008777704892428134,
          0.0013001547948194997
        ],
        [
          -0.04752564184773701,
          0.03686627308173653,
          0.0013001547948194867,
          0.006815632429983213
        ]
      ],
      "parameterCorrelation": -0.7599668610370407
    },
    "branched": {
      "parameters": [
        0.369140625,
        -0.01562402006420971,
        -0.047998046875000006,
        -0.014208984375000002
      ],
      "feederResistance": 1.446491002170128,
      "mixingLoss": 0.9844974017491751,
      "logDiscrepancy": [
        -0.047998046875000006,
        -0.014208984375000002
      ],
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "covariance": [
        [
          0.06121176769382574,
          0.0563476734679427,
          0.013709016643888962,
          0.011035228241072802
        ],
        [
          0.05634767346794267,
          0.07964592636821281,
          0.016946255330678377,
          0.013679469662845858
        ],
        [
          0.013709016643888957,
          0.01694625533067838,
          0.004014497367063168,
          0.003019937891541659
        ],
        [
          0.011035228241072798,
          0.013679469662845858,
          0.003019937891541659,
          0.0027060701131147277
        ]
      ],
      "parameterCorrelation": 0.8070059644952623
    }
  },
  "summary": {
    "feederRMS": 13.303836876205887,
    "outletRMS": 8.58930696241025,
    "coverage": 0.6818181818181818,
    "passed": false
  },
  "excluded": [
    {
      "id": "cad45-ds10mm-p3bar-sst-v1",
      "reason": "CFD acceptance failed"
    },
    {
      "id": "cad45-dm10mm-p3bar-sst-v1",
      "reason": "CFD acceptance failed"
    },
    {
      "id": "cad45-dm40mm-p3bar-sst-v1",
      "reason": "CFD acceptance failed"
    },
    {
      "id": "annular-aequal-p3bar-sst-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad90-ball-area00-sst-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad90-ball-area25-sst-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad90-ball-area50-sst-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad90-ball-area75-sst-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad45-ball00deg-p3bar-sst-mesh2-cpu64-v3-postproc-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad45-ball30deg-p3bar-sst-mesh2-cpu64-v3-postproc-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad45-ball60deg-p3bar-sst-mesh2-cpu64-v3-postproc-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad45-n01-axial-fixedmj-sst-window-v2",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad45-n04-axial-fixedmj-sst-window-v2",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad45-n04-inward15-fixedmj-sst-budget-v3",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad45-n04-inward30-fixedmj-sst-budget-v4",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "coflow-35-sst-v2",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "air-film-35-sst-v3",
      "reason": "CFD acceptance failed"
    },
    {
      "id": "early-dilution-35-sst-v1",
      "reason": "CFD acceptance failed"
    }
  ],
  "limitations": [
    "Only two accepted British nozzle cases; family uncertainty and extrapolation are weakly identified.",
    "Scalar and particle validation are separate from gas-flow cross-validation.",
    "Model discrepancy is for prediction intervals/comparisons, not a correction to conservation-balanced live flow.",
    "Single-mesh gas CFD cannot establish vacuum, powder-loaded or particle accuracy."
  ]
}
;
'use strict';
/* STEP triangle signed distance, with a lazy 0.1 mm trilinear grid (0.05 mm for the short dual-suction case). The native
 * surface supplies distance/sign; the existing inlet/outlet planes remain open.
 * No invalid STEP is auto-repaired. Fluid-volume certification is separate. */
const EductorWalls=(()=>{
 const DATA=typeof EDUCTOR_WALLS!=='undefined'?EDUCTOR_WALLS:typeof require==='function'?require('../phases/walls.json'):{};
 const cache=new Map(),sub=(a,b)=>a.map((v,i)=>v-b[i]),dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
 const bytes=s=>typeof Buffer!=='undefined'?new Uint8Array(Buffer.from(s,'base64')):Uint8Array.from(atob(s),c=>c.charCodeAt(0));
 function nearest(p,a,b,c){
  const ab=sub(b,a),ac=sub(c,a),ap=sub(p,a),d1=dot(ab,ap),d2=dot(ac,ap);if(d1<=0&&d2<=0)return a;
  const bp=sub(p,b),d3=dot(ab,bp),d4=dot(ac,bp);if(d3>=0&&d4<=d3)return b;
  const vc=d1*d4-d3*d2;if(vc<=0&&d1>=0&&d3<=0)return a.map((v,i)=>v+d1/(d1-d3)*ab[i]);
  const cp=sub(p,c),d5=dot(ab,cp),d6=dot(ac,cp);if(d6>=0&&d5<=d6)return c;
  const vb=d5*d2-d1*d6;if(vb<=0&&d2>=0&&d6<=0)return a.map((v,i)=>v+d2/(d2-d6)*ac[i]);
  const va=d3*d6-d5*d4;if(va<=0&&d4-d3>=0&&d5-d6>=0){const w=(d4-d3)/(d4-d3+d5-d6);return b.map((v,i)=>v+w*(c[i]-v));}
  const v=vb/(va+vb+vc),w=vc/(va+vb+vc);return a.map((q,i)=>q+ab[i]*v+ac[i]*w);
 }
 function build(key){
  if(cache.has(key))return cache.get(key);const data=DATA[key];if(!data)return null;
  const vb=bytes(data.v),tb=bytes(data.t),vs=new Uint16Array(vb.buffer,vb.byteOffset,vb.byteLength/2),ts=new Uint32Array(tb.buffer,tb.byteOffset,tb.byteLength/4),v=Array.from({length:vs.length/3},(_,i)=>[0,1,2].map(k=>(data.lo[k]+data.span*vs[i*3+k]/65535)/1000));
  const tris=Array.from({length:ts.length/3},(_,i)=>{const a=v[ts[i*3]],b=v[ts[i*3+1]],c=v[ts[i*3+2]],n=cross(sub(b,a),sub(c,a)),length=Math.hypot(...n);return {a,b,c,n:n.map(v=>v/(length||1)),lo:a.map((q,k)=>Math.min(q,b[k],c[k])),hi:a.map((q,k)=>Math.max(q,b[k],c[k]))};});
  function node(list){const lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity];for(const t of list)for(let k=0;k<3;k++){lo[k]=Math.min(lo[k],t.lo[k]);hi[k]=Math.max(hi[k],t.hi[k]);}if(list.length<=10)return {lo,hi,list};const span=hi.map((q,k)=>q-lo[k]),axis=span.indexOf(Math.max(...span));list.sort((a,b)=>a.lo[axis]+a.hi[axis]-b.lo[axis]-b.hi[axis]);const mid=list.length>>1;return {lo,hi,left:node(list.slice(0,mid)),right:node(list.slice(mid))};}
  const root=node(tris),grid=new Map(),h=key==='nov5'?.00005:.0001;
  const box=(p,n)=>p.reduce((s,q,i)=>s+(q<n.lo[i]?n.lo[i]-q:q>n.hi[i]?q-n.hi[i]:0)**2,0);
  function exact(p){let best=Infinity,sign=1,normal=[0,0,0];function search(n){if(box(p,n)>best)return;if(n.list){for(const t of n.list){const q=nearest(p,t.a,t.b,t.c),d=sub(p,q),r=dot(d,d);if(r<best){best=r;sign=dot(d,t.n)>=0?1:-1;normal=t.n;}}}else{const l=box(p,n.left),r=box(p,n.right);if(l<r){search(n.left);search(n.right);}else{search(n.right);search(n.left);}}}search(root);
   // Oriented ray winding gives the sign of the UNION of closed hardware
   // solids. A nearest-face normal alone is discontinuous at concave edges
   // and at assembly-part overlaps; it is not a valid signed-distance sign.
   const direction=[.9238795325,.331413574,.191341716];let winding=0;
   function raybox(n){let lo=0,hi=Infinity;for(let i=0;i<3;i++){const a=(n.lo[i]-p[i])/direction[i],b=(n.hi[i]-p[i])/direction[i];lo=Math.max(lo,Math.min(a,b));hi=Math.min(hi,Math.max(a,b));}return hi>=lo;}
   function ray(n){if(!raybox(n))return;if(n.list){for(const t of n.list){const e1=sub(t.b,t.a),e2=sub(t.c,t.a),h=cross(direction,e2),det=dot(e1,h);if(Math.abs(det)<1e-18)continue;const q=sub(p,t.a),u=dot(q,h)/det;if(u<0||u>1)continue;const v=cross(q,e1),w=dot(direction,v)/det;if(w<0||u+w>1)continue;const distance=dot(e2,v)/det;if(distance>1e-10)winding+=dot(direction,t.n)>0?1:-1;}}else{ray(n.left);ray(n.right);}}
   ray(root);sign=winding!==0?-1:1;return {distance:sign*Math.sqrt(best),normal};}
  function distance(x,y,z){const q=[x/h,y/h,z/h],i=q.map(Math.floor),f=q.map((v,k)=>v-i[k]);let value=0,reference=null,sharp=false;for(let mask=0;mask<8;mask++){const j=i.map((v,k)=>v+((mask>>k)&1)),key=j.join(','),weight=f.reduce((s,v,k)=>s*(((mask>>k)&1)?v:1-v),1);let d=grid.get(key);if(d===undefined){d=exact(j.map(v=>v*h));if(grid.size>=100000)grid.delete(grid.keys().next().value);grid.set(key,d);}if(reference&&dot(reference,d.normal)<.99)sharp=true;reference??=d.normal;value+=weight*d.distance;}return sharp?exact([x,y,z]).distance:value;}
  function normal(x,y,z){const q=[x,y,z],d=q.map((_,i)=>{const a=q.slice(),b=q.slice();a[i]-=h/2;b[i]+=h/2;return distance(...b)-distance(...a);}),n=Math.hypot(...d);return d.map(v=>v/(n||1));}
  const out={distance,exact,normal,gridSpacing:h,source:data,triangles:tris.length};if(cache.size>=2)cache.delete(cache.keys().next().value);cache.set(key,out);return out;
 }
 function attach(geo,key){const sdf=build(key);if(!sdf)return {...geo,wallStatus:'No valid native STEP SDF; reduced walls retained.'};
  const clearance=(x,y,dp=0,z=0)=>{const inlet=Math.max(...geo.branches.map(b=>{const a=b.segments[0];return dot(sub([x,y,z],a.a),a.tangent);}));return Math.min(sdf.distance(x,y,z),geo.end-x,inlet)-dp/2;};
  return {...geo,clearance,sdf,wallStatus:'Native STEP surface distance · '+(sdf.gridSpacing*1000)+' mm grid · fluid-volume gate pending',wallNormal:sdf.normal};
 }
 return {build,attach,nearest,register:(key,data)=>{DATA[key]=data;cache.delete(key);},has:key=>!!DATA[key]};
})();
if(typeof module!=='undefined')module.exports=EductorWalls;

'use strict';
// Predictor/corrector pseudo-arclength continuation in scaled (P, boundary P).
const EductorContinuation=(()=>{
 function trace(residual,a,b,{step=.005,count=80,tolerance=1e-9}={}){
  const points=[a.slice(),b.slice()];let ds=step;
  for(let n=0;n<count;n++){
   const u=points.at(-2),v=points.at(-1),delta=v.map((x,i)=>x-u[i]),norm=Math.hypot(...delta);if(!norm)break;
   const t=delta.map(x=>x/norm),pred=v.map((x,i)=>x+ds*t[i]);let x=pred.slice(),good=false;
   for(let k=0;k<16;k++){
    const f=residual(...x),arc=(x[0]-pred[0])*t[0]+(x[1]-pred[1])*t[1];if(!Number.isFinite(f))break;
    if(Math.abs(f)<tolerance&&Math.abs(arc)<tolerance){good=true;break;}
    const h=1e-6,fx=(residual(x[0]+h,x[1])-residual(x[0]-h,x[1]))/(2*h),fy=(residual(x[0],x[1]+h)-residual(x[0],x[1]-h))/(2*h),det=fx*t[1]-fy*t[0];if(!Number.isFinite(det)||Math.abs(det)<1e-12)break;
    x=[x[0]-(f*t[1]-fy*arc)/det,x[1]-(fx*arc-f*t[0])/det];
   }
   if(good){points.push(x);ds=Math.min(step,ds*1.4);}else{ds*=.5;if(ds<step/128)break;n--;}
  }
  return points;
 }
 return {trace};
})();
if(typeof module!=='undefined')module.exports=EductorContinuation;

'use strict';
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
const ExposurePopulation=(()=>{
 const sum=a=>a.reduce((s,v)=>s+v,0);
 function create(diameters,shares,density,{bins=36,ratio=2}={}){
  const mass=d=>density*Math.PI*d**3/6,min=Math.min(...diameters.map(mass)),max=Math.max(...diameters.map(mass))*128,r=Math.max(ratio,(max/min)**(1/(bins-1))),m=Array.from({length:bins},(_,i)=>min*r**i),q=new Float64Array(bins);
  diameters.forEach((d,i)=>deposit(q,m,mass(d),shares[i]/mass(d)));
  return {m,q,density,initialNumber:sum(q),initialMass:1,overflowNumber:0,overflowMass:0,steps:0,maxMassError:0};
 }
 function deposit(q,m,v,n){if(v>m.at(-1))return {number:n,mass:n*v};let k=0;while(k<m.length-2&&m[k+1]<v)k++;const w=(m[k+1]-v)/(m[k+1]-m[k]);q[k]+=n*w;q[k+1]+=n*(1-w);return {number:0,mass:0};}
 function advance(s,C,dt,kernel,alpha=1){
  if(!(C>0&&dt>0&&alpha>0))return s;
  const n=s.m.length,K=Array.from({length:n},()=>new Map()),getK=(i,j)=>{if(!K[i].has(j)){const v=Math.max(0,alpha*kernel(i,j));K[i].set(j,v);K[j].set(i,v);}return K[i].get(j);};
  const births=Array.from({length:n},()=>new Map());
  function birth(i,j){if(births[i].has(j))return births[i].get(j);const v=s.m[i]+s.m[j];let out;
   if(v>s.m.at(-1))out={bins:[],mass:v};else{let l=0;while(l<n-2&&s.m[l+1]<v)l++;const w=(s.m[l+1]-v)/(s.m[l+1]-s.m[l]);out={bins:[[l,w],[l+1,1-w]],mass:0};}
   out.lossI=(i===j?2:1)-(out.bins.find(b=>b[0]===i)?.[1]||0);out.lossJ=i===j?out.lossI:1-(out.bins.find(b=>b[0]===j)?.[1]||0);births[i].set(j,out);return out;
  }
  function react(i,j,h){const k=C*getK(i,j),a=s.q[i],b=s.q[j];if(!(a>0&&b>0&&k>0))return;const out=birth(i,j),si=out.lossI,sj=out.lossJ;let event;
   // Fixed-pivot birth may return most of a newly merged particle to its
   // parent's bin. Integrate NET stoichiometry, otherwise giant/small pairs
   // acquire an artificial collision-timescale stiffness and biased depletion.
   if(i===j){const z=.5*k*si*a*h;event=a/si*(z/(1+z));}
   else if(sj<1e-15)event=a/si*(-Math.expm1(-si*k*b*h));
   else {const aa=a/si,bb=b/sj,small=Math.min(aa,bb),large=Math.max(aa,bb),delta=large-small,kk=k*si*sj,z=kk*delta*h;event=delta<large*1e-10?small*(kk*large*h/(1+kk*large*h)):small*large*(-Math.expm1(-z))/(delta+small*(-Math.expm1(-z)));}
   s.q[i]-=si*event;if(j!==i)s.q[j]-=sj*event;
   for(const [bin,w] of out.bins)if(bin!==i&&bin!==j)s.q[bin]+=w*event;
   if(out.mass){s.overflowNumber+=event;s.overflowMass+=event*out.mass;}
  }
  // Exact positive pair reactions with symmetric (Strang) splitting. A bound
  // on bulk and materially populated bin loss rates resolves evolution without allowing negligible
  // high-mass tails to force millions of explicit stability substeps.
  let left=dt;while(left>dt*1e-12){
   // Retain sub-threshold bins in the mass ledger but do not let numerically
   // negligible tails dominate pair work. Their total mass is bounded by
   // bins × 1e-14 of initial mass; the threshold is below the mass-test budget.
   const active=[];for(let i=0;i<n;i++)if(s.q[i]*s.m[i]>1e-14*s.initialMass)active.push(i);
   let rate=0,N=sum(s.q),resolvedLoss=0;
   for(const i of active){let loss=0;for(const j of active){const pair=i<=j?birth(i,j):birth(j,i),encounter=C*getK(i,j)*s.q[j];rate+=s.q[i]*encounter;loss+=encounter*(i===j?.5*pair.lossI:i<j?pair.lossI:pair.lossJ);}if(s.q[i]*s.m[i]>1e-8*s.initialMass)resolvedLoss=Math.max(resolvedLoss,loss);}
   const h=Math.min(left,rate>0?.02*N/rate:left,resolvedLoss>0?.05/resolvedLoss:left),pairs=[];
   for(let a=0;a<active.length;a++)for(let b=a;b<active.length;b++)pairs.push([active[a],active[b]]);
   for(const [i,j] of pairs)react(i,j,h/2);for(let k=pairs.length-1;k>=0;k--)react(...pairs[k],h/2);
   s.steps++;left-=h;if(s.steps>1000000)throw Error('Population balance exceeded its substep budget');
  }
  const total=s.q.reduce((v,q,i)=>v+q*s.m[i],s.overflowMass);s.maxMassError=Math.max(s.maxMassError,Math.abs(total-s.initialMass)/s.initialMass);return s;
 }
 function summary(s){const number=sum(s.q)+s.overflowNumber,mass=s.q.reduce((v,q,i)=>v+q*s.m[i],s.overflowMass);return {numberPerKg:number,numberFraction:number/s.initialNumber,numberLoss:1-number/s.initialNumber,mass,massError:s.maxMassError,meanDiameter:Math.cbrt(6*mass/(s.density*Math.PI*number)),overflowMassFraction:s.overflowMass/mass,distribution:s.m.map((m,i)=>({diameter:Math.cbrt(6*m/(s.density*Math.PI)),numberPerKg:s.q[i],massFraction:m*s.q[i]})),steps:s.steps};}
 return {create,advance,summary};
})();
if(typeof module!=='undefined')module.exports=ExposurePopulation;

'use strict';
const EductorJet=(()=>{
 const R=287.05,cp=1004.675;
 function notional(flow,dn){
  const thrust=flow.network?.motiveThrust??flow.mp*flow.uJet,u=flow.mp?thrust/flow.mp:0,T=flow.T0-u*u/(2*cp),rho=flow.portP/(R*T),A=u>0?flow.mp/(rho*u):0;
  return {u,T,rho,A,diameter:Math.sqrt(4*A/Math.PI),thrust,power:.5*flow.mp*u*u,machDisk:dn*.67*Math.sqrt(flow.p0/flow.portP),valid:T>0&&Number.isFinite(A)&&A>0,source:'Mass + momentum + adiabatic energy equivalent source; Birch-type momentum formulation (1987). No shock-cell resolution.'};
 }
 function annulusMass(amp,a){
  const base=a.mean-amp*a.fb;if(base>=0)return 0;
  const y0=Math.min(1,-a.b*a.b*Math.log(-base/amp)/(a.r*a.r));
  return -a.rho*Math.PI*a.r*a.r*(base*(1-y0)+amp*a.b*a.b/(a.r*a.r)*(Math.exp(-y0*a.r*a.r/(a.b*a.b))-Math.exp(-a.r*a.r/(a.b*a.b))));
 }
 // One gas-only inlet diagnostic per passage: never recompute onset from
 // the expanding local area. Powder affects the solved flow, not this gas Ct.
 function confinement(f,r){
  const rho=f.portP/(R*f.Tin),rhoA=rho*Math.PI*r*r,uk=f.mg/rhoA;
  const momentum=(f.network?.motiveThrust??f.mp*f.uJet)+f.ms*(f.network?.secondaryVelocity||0)+f.mh*(f.network?.bleedVelocity||0),ud2=momentum/rhoA;
  return {Ct:ud2>.5*uk*uk?uk/Math.sqrt(ud2-.5*uk*uk):null,rho,uk,ud2,radius:r,convention:'Gas-only, fixed reference bore; all supplied gas mass and axial momentum. Candidate threshold, not universal.'};
 }
 function recirculation(a,f,dEq,spread,reference){
  const fixed=reference||confinement(f,a.r),Ct=fixed.Ct;
  const rhoA=a.rho*Math.PI*a.r*a.r,ud2=((f.network?.motiveThrust??f.mp*f.uJet)+f.ms*(f.network?.secondaryVelocity||0)+f.mh*(f.network?.bleedVelocity||0))/rhoA;
  // The cylindrical confined-jet closure has no validated diffuser/junction
  // extension. Keep its diagnostic, but do not impose reverse flow there.
  const supported=fixed.constantArea!==false;
  const entrained=.32*f.mp*a.x/Math.max(dEq,1e-9),supply=f.ms*a.F+f.mh*a.H,closure=Math.max(0,Math.min(1,(1-a.b/a.r)/.2)),demand=supported&&Ct!==null&&Ct<.98?Math.max(0,entrained-supply)*closure:0;
  if(!demand)return {amp:a.amp,Ct,supported,reason:supported?null:'Variable-area passage: diffuser recirculation not modeled',demand:0,reverse:annulusMass(a.amp,a),momentumLimited:false};
  let lo=a.mean/Math.max(a.fb,1e-12),hi=Math.max(lo*2,a.amp);while(annulusMass(hi,a)<demand&&hi<1e6)hi*=2;
  for(let k=0;k<45;k++){const mid=(lo+hi)/2;if(annulusMass(mid,a)<demand)lo=mid;else hi=mid;}
  const variance=a.b*a.b/(2*a.r*a.r)*(-Math.expm1(-2*a.r*a.r/(a.b*a.b)))-a.fb*a.fb;
  const budget=Math.max(0,ud2-a.mean*a.mean),cap=Math.sqrt(budget/Math.max(variance,1e-16)),amp=Math.max(a.amp,Math.min((lo+hi)/2,cap));
  return {amp,Ct,supported,demand,reverse:annulusMass(amp,a),momentumLimited:cap<(lo+hi)/2,closureX:(a.r-dEq/Math.sqrt(2))/Math.max(spread,1e-9)};
 }
 return {notional,annulusMass,confinement,recirculation};
})();
if(typeof module!=='undefined')module.exports=EductorJet;

const EDUCTOR_TOPOLOGY={"nov0":{"source":"British Eductor 1'' - w-o converge-diverge.step","sha256":"b9c8c52bc3705774c7eacd6cc4ef87c254e0533a7784cb32dbcc6c931af8f95a","coordinates":"nozzle-relative CAD mm; flow +x","cadValid":true,"common":{"points":[[11.301700000000004,-8.938606635269468e-16,38.5],[11.301700000000004,-8.938606635269468e-16,38.5]],"radius":8.0,"lengthMM":0.0},"branches":[{"face":11,"surface":"CYLINDER","radius":8.0,"points":[[11.301700000000004,-8.938606635269468e-16,38.5],[11.301700000000004,0.0,0.0]],"maxCircleResidualMM":0,"junctionContinuationMM":0,"id":0,"lengthMM":38.5,"turnRadians":0.0,"nativeCenterChecks":0,"centersInsideSolid":0}],"bleedCount":8,"notes":"Internal swept/cylindrical centreline extraction. Merger and main-bore continuations are reduced tube unions, not a CAD-resolved flow field."},"nov1":{"source":"British Eductor 1'' - with converge-diverge.step","sha256":"37591a07f7a1e3917abab0396dc8b8592e9b8c2540d445da05dd5a66dd0b497c","coordinates":"nozzle-relative CAD mm; flow +x","cadValid":true,"common":{"points":[[11.301700000000004,-8.938606635269468e-16,38.5],[11.301700000000004,-8.938606635269468e-16,38.5]],"radius":8.0,"lengthMM":0.0},"branches":[{"face":32,"surface":"CYLINDER","radius":8.0,"points":[[11.301700000000004,-8.938606635269468e-16,38.5],[11.301700000000004,0.0,0.0]],"maxCircleResidualMM":0,"junctionContinuationMM":0,"id":0,"lengthMM":38.5,"turnRadians":0.0,"nativeCenterChecks":0,"centersInsideSolid":0}],"bleedCount":8,"notes":"Internal swept/cylindrical centreline extraction. Merger and main-bore continuations are reduced tube unions, not a CAD-resolved flow field."},"nov2":{"source":"British Eductor - original size - with converge-diverge.step","sha256":"f4c5f9e89805e8536577b542324c9ef3095e631ca88cc4650f3c3275c808db77","coordinates":"nozzle-relative CAD mm; flow +x","cadValid":true,"common":{"points":[[1.8560000000000016,-8.938606635269468e-16,38.5],[1.8560000000000016,-8.938606635269468e-16,38.5]],"radius":5.0,"lengthMM":0.0},"branches":[{"face":32,"surface":"CYLINDER","radius":5.0,"points":[[1.8560000000000016,-8.938606635269468e-16,38.5],[1.8560000000000016,0.0,0.0]],"maxCircleResidualMM":0,"junctionContinuationMM":0,"id":0,"lengthMM":38.5,"turnRadians":0.0,"nativeCenterChecks":0,"centersInsideSolid":0}],"bleedCount":8,"notes":"Internal swept/cylindrical centreline extraction. Merger and main-bore continuations are reduced tube unions, not a CAD-resolved flow field."},"nov3":{"source":"British Eductor - original size w-o converge-diverge.step","sha256":"4a72e84eec0c32029692e3d7335a547c0d571fd3f80a0fcca6915538b0fa36fe","coordinates":"nozzle-relative CAD mm; flow +x","cadValid":true,"common":{"points":[[1.8560000000000016,-8.938606635269468e-16,38.5],[1.8560000000000016,-8.938606635269468e-16,38.5]],"radius":5.0,"lengthMM":0.0},"branches":[{"face":57,"surface":"CYLINDER","radius":5.0,"points":[[1.8560000000000016,-8.938606635269468e-16,38.5],[1.8560000000000016,0.0,0.0]],"maxCircleResidualMM":0,"junctionContinuationMM":0,"id":0,"lengthMM":38.5,"turnRadians":0.0,"nativeCenterChecks":0,"centersInsideSolid":0}],"bleedCount":8,"notes":"Internal swept/cylindrical centreline extraction. Merger and main-bore continuations are reduced tube unions, not a CAD-resolved flow field."},"nov4":{"source":"0.75'' - SHORT - 3D printed.step","sha256":"831ec6fd67f0e776d9952f5b0efdd555a9312f04215daa8fb6acbb2ce2046b55","coordinates":"nozzle-relative CAD mm; flow +x","cadValid":false,"common":{"points":[[-27.782294321497375,-8.93860663526947e-16,50.34600282048214],[-27.782294321497375,-8.93860663526947e-16,50.34600282048214]],"radius":8.0,"lengthMM":0.0},"branches":[{"face":10,"surface":"CYLINDER","radius":8.0,"points":[[-27.782294321497375,-8.93860663526947e-16,50.34600282048214],[22.56370849898483,0.0,0.0]],"maxCircleResidualMM":0,"junctionContinuationMM":0,"id":0,"lengthMM":71.19999999999999,"turnRadians":0.0,"nativeCenterChecks":0,"centersInsideSolid":null}],"bleedCount":0,"notes":"Internal swept/cylindrical centreline extraction. Merger and main-bore continuations are reduced tube unions, not a CAD-resolved flow field. Printed tube source is invalid and requires CAD repair before CFD."},"nov5":{"source":"1'' Eductor dual suction SHORT.step","sha256":"2f5a509a8bb37b38b9f4f648afff39246f9455002406ad875de0ac547ca52d07","coordinates":"nozzle-relative CAD mm; flow +x","cadValid":true,"common":{"points":[[-91.08927522873023,-65.56014756989518,2.60140027202554e-15],[-62.610727125275154,-37.35946474537766,2.0898733437978604e-15]],"radius":7.95,"lengthMM":40.07875015328999},"branches":[{"face":112,"surface":"TORUS","radius":7.950000000000001,"points":[[-62.610727125275154,-37.35946474537766,2.0898733437978604e-15],[-45.15613312638333,-20.07517527228627,27.103175552614033],[-44.65115546415441,-19.582505808291558,27.593245500883302],[-44.13157426762722,-19.09013310995987,28.06811253538099],[-43.5977025130862,-18.598353764193657,28.527490613961945],[-43.049861784913524,-18.107464000482025,28.971103024448837],[-42.48838208187876,-17.617759512463067,29.39868255131337],[-41.91360161835991,-17.129535279809133,29.809971636636785],[-41.32586662061578,-16.643085390542076,30.204722535253175],[-40.72553111823219,-16.15870286388568,30.58269746398191],[-40.1129567308678,-15.676679473761771,30.943668744859274],[-39.48851245042794,-15.197305573036553,31.287418942283264],[-38.852574418797595,-14.720869918622837,31.61374099398864],[-38.20552570126785,-14.247659497543829,31.922438335773396],[-37.5477560557914,-13.777959354062663,32.21332501990175],[-36.87966169820708,-13.312052417982558,32.48622582711219],[-36.20164506357427,-12.850219334220448,32.74097637216284],[-35.51411456376095,-12.392738293756908,32.977423202851085],[-34.81748434143183,-11.93988486606443,33.19542389244748],[-34.112174020584156,-11.49193183311456,33.39484712548814],[-33.398608453781954,-11.049149025064297,33.5755727768744],[-32.677217466240684,-10.611803157720482,33.73749198423165],[-31.948435596916745,-10.180157671880236,33.88050721348395],[-31.212701836757283,-9.754472574643907,34.004532317604834],[-30.470459364268557,-9.335004282796698,34.1094925885092],[-29.72215527856182,-8.922005468352747,34.19532480205455],[-28.968240330037354,-8.515724906354768,34.26197725612496],[-28.20916864886943,-8.116407325021532,34.309409801774365],[-27.445397471455358,-7.724293258332619,34.33759386741088],[-26.677386864993075,-7.339618901139672,34.346512476007234],[-25.905599450354195,-6.962615966891581,34.33616025532713],[-25.130500123418308,-6.593511548058892,34.30654344116117],[-24.352555775037096,-6.2325279793418,34.257679873570794],[-23.57223500979685,-5.879882703744053,34.18959898614199],[-22.790007863748436,-5.535788141593412,34.10234178825561],[-22.00634552127527,-5.200451562587587,33.9959608403849],[-21.22172003126952,-4.874074960942755,33.870520222435076],[-20.4366040227873,-4.55685493371975,33.72609549514381],[-19.651470420354663,-4.248982562401294,33.56277365456645],[-18.866792159095745,-3.950643297791816,33.38065307967276],[-18.083041899854088,-3.6620168483086886,33.179843473087125],[-17.30069174447963,-3.383277071732657,32.96046579500783],[-16.52021295145186,-3.114591870482349,32.72265219034504],[-15.742075652011827,-2.8561230904762107,32.46654590912179],[-14.966748566972278,-2.608026423642576,32.19230122018549],[-14.194698724377893,-2.3704513141366204,31.900083318282096],[-13.426391178185021,-2.1435408683207684,31.590068224548933],[-12.662288728130626,-1.927431768562857,31.262442680486263],[-11.902851640959007,-1.7222541909037616,30.917404035471154],[-11.148537373174243,-1.528131726644162,30.555160127881496],[-10.399800295485697,-1.3451813078978385,30.175929159902097],[-9.657091419111829,-1.1735131371560663,29.779939566087656],[-8.920858124108143,-1.0132306209058466,29.367429875762557],[-8.191543889881949,-0.8644303073415854,28.938648569339435],[-7.469588028057167,-0.7272018282081875,28.49385392864417],[-6.7554254178493025,-0.6016278448101776,28.033313881336355],[-6.049486244110469,-0.487783998219643,27.557305839519803],[-5.352195738202141,-0.3857388637128596,27.066116532639754],[-4.66397392185177,-0.29555390946316257,26.560041834768036],[-3.9852353541471714,-0.21728345951474456,26.039386586379166],[-3.316388881821922,-0.15097466105994314,25.504464410726314],[-2.657837392981122,-0.09666745603949249,24.955597524926002],[-2.009977574416858,-0.05439455708302911,24.393116545866864],[-1.3731996726588278,-0.024181427804192825,23.817360291058097],[-0.7478872589042922,-0.006046267462319885,23.228675574538677],[-0.1344169979694385,1.126553516914366e-13,22.627416997969537],[22.950233055063173,0.22751871995337655,0.0022879324922016053]],"maxCircleResidualMM":1.543649238600926e-14,"junctionContinuationMM":32.324128728682155,"id":0,"lengthMM":123.8791882809795,"turnRadians":1.7736926163191722,"nativeCenterChecks":16,"centersInsideSolid":0},{"face":115,"surface":"TORUS","radius":7.950000000000001,"points":[[-62.610727125275154,-37.35946474537766,2.0898733437978604e-15],[-45.15613312638333,-20.075175272286277,-27.10317555261403],[-44.65115546415442,-19.582505808291593,-27.593245500883288],[-44.13157426762724,-19.09013310995988,-28.068112535380983],[-43.59770251308624,-18.598353764193693,-28.527490613961913],[-43.049861784913524,-18.107464000482025,-28.971103024448837],[-42.488382081878775,-17.61775951246309,-29.398682551313357],[-41.91360161835991,-17.129535279809133,-29.809971636636785],[-41.325866620615805,-16.64308539054209,-30.204722535253165],[-40.72553111823219,-16.158702863885676,-30.58269746398191],[-40.112956730867815,-15.676679473761782,-30.943668744859256],[-39.48851245042794,-15.19730557303655,-31.287418942283267],[-38.85257441879761,-14.720869918622842,-31.61374099398864],[-38.205525701267895,-14.247659497543856,-31.92243833577337],[-37.5477560557914,-13.777959354062665,-32.213325019901745],[-36.87966169820711,-13.31205241798258,-32.48622582711218],[-36.20164506357427,-12.850219334220451,-32.74097637216284],[-35.51411456376098,-12.392738293756924,-32.97742320285108],[-34.81748434143183,-11.93988486606443,-33.19542389244748],[-34.11217402058417,-11.491931833114561,-33.394847125488134],[-33.39860845378194,-11.04914902506429,-33.57557277687442],[-32.6772174662407,-10.611803157720496,-33.73749198423166],[-31.948435596916802,-10.180157671880261,-33.880507213483945],[-31.212701836757283,-9.75447257464391,-34.00453231760483],[-30.4704593642686,-9.335004282796735,-34.109492588509184],[-29.722155278561836,-8.922005468352765,-34.19532480205454],[-28.968240330037382,-8.515724906354786,-34.261977256124965],[-28.20916864886943,-8.116407325021532,-34.30940980177437],[-27.4453974714554,-7.724293258332646,-34.337593867410874],[-26.677386864993103,-7.339618901139677,-34.346512476007234],[-25.905599450354195,-6.962615966891583,-34.33616025532713],[-25.130500123418294,-6.593511548058885,-34.306543441161175],[-24.35255577503711,-6.232527979341807,-34.257679873570794],[-23.572235009796835,-5.879882703744045,-34.18959898614199],[-22.79000786374845,-5.535788141593421,-34.1023417882556],[-22.006345521275314,-5.200451562587604,-33.99596084038492],[-21.221720031269534,-4.874074960942766,-33.87052022243508],[-20.436604022787314,-4.5568549337197535,-33.72609549514382],[-19.651470420354663,-4.248982562401294,-33.56277365456645],[-18.86679215909576,-3.950643297791825,-33.38065307967276],[-18.083041899854088,-3.662016848308689,-33.179843473087125],[-17.30069174447963,-3.3832770717326572,-32.96046579500783],[-16.52021295145186,-3.1145918704823417,-32.72265219034503],[-15.742075652011856,-2.8561230904762107,-32.46654590912179],[-14.96674856697232,-2.608026423642589,-32.19230122018551],[-14.194698724377908,-2.3704513141366266,-31.900083318282103],[-13.426391178185064,-2.143540868320771,-31.59006822454894],[-12.662288728130626,-1.9274317685628566,-31.262442680486263],[-11.902851640959021,-1.7222541909037723,-30.91740403547117],[-11.148537373174243,-1.5281317266441607,-30.555160127881496],[-10.399800295485726,-1.3451813078978514,-30.17592915990212],[-9.657091419111829,-1.1735131371560654,-29.779939566087656],[-8.920858124108143,-1.0132306209058515,-29.367429875762568],[-8.191543889882006,-0.8644303073416015,-28.938648569339488],[-7.469588028057181,-0.7272018282081895,-28.49385392864418],[-6.755425417849274,-0.6016278448101723,-28.033313881336337],[-6.049486244110476,-0.48778399821964763,-27.557305839519813],[-5.352195738202184,-0.38573886371286126,-27.06611653263978],[-4.66397392185177,-0.29555390946316257,-26.560041834768036],[-3.9852353541472,-0.21728345951474853,-26.039386586379187],[-3.316388881821922,-0.15097466105994314,-25.504464410726314],[-2.6578373929811505,-0.09666745603949538,-24.955597524926034],[-2.0099775744168724,-0.05439455708302548,-24.393116545866864],[-1.373199672658835,-0.024181427804191472,-23.817360291058105],[-0.7478872589042851,-0.006046267462318323,-23.228675574538663],[-0.13441699796944562,1.2422502172768851e-13,-22.627416997969554],[22.95023305506387,0.2275187199537769,-0.0022879324922122635]],"maxCircleResidualMM":1.56836440602695e-14,"junctionContinuationMM":32.324128728682666,"id":1,"lengthMM":123.87918828098002,"turnRadians":1.7736926163192184,"nativeCenterChecks":16,"centersInsideSolid":0}],"bleedCount":0,"notes":"Internal swept/cylindrical centreline extraction. Merger and main-bore continuations are reduced tube unions, not a CAD-resolved flow field."},"nov6":{"source":"1'' Eductor 4 suction.step","sha256":"cc059c87f555a10f6d866e7317e8a388f83d9d7866b57793e737025241618f60","coordinates":"nozzle-relative CAD mm; flow +x","cadValid":true,"common":{"points":[[-114.70979579551081,-85.30330085889939,-2.4477358023231794e-16],[-95.76202633310385,-65.89376104429411,-4.98168173379554e-12]],"radius":7.0,"lengthMM":27.1246788628966},"branches":[{"face":98,"surface":"BSPLINE","radius":7.000063570589805,"points":[[-95.76202633310385,-65.89376104429411,-4.98168173379554e-12],[-98.19731008059662,-66.38946973937996,1.7183078868366544],[-97.31773537238979,-64.93034520939152,2.1389982829491254],[-96.4152482480099,-63.39833421705697,2.5982253214630395],[-95.49000282895807,-61.7972274114214,3.0938208869478605],[-94.54215323676488,-60.13081544171301,3.6236168638880053],[-93.57185359296088,-58.40288895715984,4.185445136768018],[-92.57925801907655,-56.61723860699022,4.777137590072307],[-91.56452063664241,-54.777655040432144,5.39652610828541],[-90.52779556718906,-52.88792890671376,6.041442575891776],[-89.46923693224701,-50.95185085506331,6.70971887737593],[-88.38899885334669,-48.97321153470846,7.399186897222471],[-87.28723545200592,-46.95580159477783,8.107678519959224],[-86.16410084974163,-44.903411684417364,8.833025630106558],[-85.01974916807873,-42.81983245284562,9.573060112153374],[-83.85433452854227,-40.7088545492809,10.325613850588704],[-82.66801105265708,-38.57426862294178,11.088518729901416],[-81.46093286194824,-36.41986532304659,11.85960663458053],[-80.23325407794061,-34.24943529881389,12.636709449114932],[-78.98512882215928,-32.066769199461945,13.417659057993635],[-77.71671121612911,-29.87565767420878,14.200287345705824],[-76.42815538139695,-27.679891372138194,14.982426196795714],[-75.11961543952552,-25.48326094250075,15.761907495731526],[-73.79124551199624,-23.28955703454075,16.536563127000996],[-72.44319972030158,-21.102570297437563,17.304224975105413],[-71.07563218603022,-18.926091380435157,18.06272492452032],[-69.68869703068214,-16.76391093266053,18.809894859785317],[-68.28254837577211,-14.619819603315245,19.543566665393868],[-66.85734034281609,-12.49760804164692,20.26157222581189],[-65.41322705333002,-10.40106689690299,20.96174342550539],[-63.95036262886541,-8.333986818246053,21.64191214900043],[-62.468901190971714,-6.300158454831156,22.299910280826797],[-60.968996861170716,-4.303372455873843,22.93356970547074],[-59.45080376098418,-2.3474194705898066,23.540722307418367],[-57.91447601193389,-0.43609014819460373,24.11919997115588],[-56.36016773554167,1.4268248620962192,24.666834581169446],[-54.78803305332932,3.2375349110670646,25.18145802194528],[-53.1982260868186,4.992249349502285,25.66090217796951],[-51.59090095753131,6.687177528186497,26.10299893372843],[-49.96621178699163,8.318528798013105,26.505580173779897],[-48.32431269673567,9.882512509817218,26.866477782647717],[-46.665357808302964,11.375338014319336,27.183523644782326],[-44.989501243233015,12.793214662240329,27.454549644634287],[-43.29689712305295,14.132351804354357,27.67738766668317],[-41.58769956923653,15.388958791660972,27.849869595532443],[-39.86206270331433,16.559244974908157,27.96982731565124],[-38.120140646826734,17.63941970480128,28.035092711485667],[-36.36208752131409,18.625692332045904,28.043497667481862],[-34.58805744831039,19.514272207570524,27.992874068211126],[-32.798204549335004,20.301368682279133,27.881053798227814],[-30.992682945905877,20.98319110687126,27.705868741970107],[-29.171646759540884,21.555948832046465,27.465150783876226],[-27.33525011176328,22.015851208562985,27.156731808418943],[-25.48364712412208,22.359107587449174,26.778443700230937],[-23.616991918143945,22.581927319458412,26.328118343782897],[-21.735438615351022,22.680519755293545,25.803587623515845],[-19.83914133726544,22.65109424565738,25.202683423870816],[-17.928254205399746,22.48986014145674,24.52323762940321],[-16.002931341280025,22.193026793599593,23.76308212467055],[-14.06332686644803,21.756803552804644,22.920048794125382],[-12.109594902445416,21.17739976979065,21.991969522220103],[-10.14188957080642,20.451024795305916,20.976676193423526],[-8.160364993021012,19.57388798027185,19.87200069229928],[-6.165175290605028,18.542198675511095,18.675774903356757],[-4.156474585080517,17.35216623182127,17.385830711091707],[-2.1344169979695238,15.999999999999995,16.000000000000007],[21.494558401330476,0.1991134025426362,-0.1942765570616416]],"maxCircleResidualMM":0.0007088845430021025,"junctionContinuationMM":32.71469224817847,"id":0,"lengthMM":184.42095652457203,"turnRadians":3.85624180686794,"nativeCenterChecks":16,"centersInsideSolid":0},{"face":99,"surface":"BSPLINE","radius":7.000282028264419,"points":[[-95.76202633310385,-65.89376104429411,-4.98168173379554e-12],[-93.32674258561471,-65.39805234921282,-1.4599234622906825],[-92.0042138408401,-64.40691896419946,-1.779486517280696],[-90.66693906386209,-63.43097846399638,-2.1253356049288743],[-89.31547455064768,-62.46983255819278,-2.4960569243234487],[-87.95037659716377,-61.523082956377316,-2.890236674552277],[-86.57220149937757,-60.59033136813875,-3.306461054703177],[-85.18150555325613,-59.67117950306581,-3.743316263863973],[-83.7788450547666,-58.765229070747296,-4.199388501122454],[-82.36477629987601,-57.87208178077186,-4.673263965566484],[-80.93985558455157,-56.99133934272834,-5.163528856283831],[-79.50463920474809,-56.12260346621677,-5.668769372382612],[-78.05968345640011,-55.265475860857215,-6.187571713006479],[-76.60554463546978,-54.4195582362445,-6.718522077254426],[-75.14277903792072,-53.58445230197206,-7.260206664223095],[-73.67194295971657,-52.75975976763336,-7.811211673009063],[-72.19359269682113,-51.94508234282193,-8.370123302708858],[-70.70828454519781,-51.14002173713114,-8.935527752419182],[-69.2165748008104,-50.3441796601545,-9.506011221236557],[-67.7190197596226,-49.5571578214855,-10.080159908257553],[-66.21617571759796,-48.778557930717554,-10.656560012578815],[-64.7085989707002,-48.00798169744414,-11.23379773329691],[-63.19684581489293,-47.24503083125876,-11.810459269508419],[-61.68147254613984,-46.489307041754856,-12.385130820309932],[-60.163035460404565,-45.740412038525875,-12.956398584798041],[-58.64209085365075,-44.99794753116528,-13.522848762069337],[-57.11919502184204,-44.26151522926653,-14.083067551220436],[-55.594904260942194,-43.530716842423146,-14.635641151347855],[-54.06977486691475,-42.80515408022855,-15.179155761548255],[-52.54436313572251,-42.08442865227853,-15.712197580923467],[-51.01922536332308,-41.36814226818104,-16.233352808602717],[-49.49491784567838,-40.655896637532315,-16.741207643688114],[-47.971996878751334,-39.94729346992616,-17.2343482852761],[-46.451018758504695,-39.24193447495634,-17.711360932463137],[-44.9325397809014,-38.539421362216665,-18.17083178434571],[-43.41711624190416,-37.83935584130087,-18.61134704002032],[-41.90530443747591,-37.14133962180277,-19.03149289858342],[-40.397660663579394,-36.444974413316125,-19.429855559131518],[-38.89474121617757,-35.74986192543475,-19.805021220761017],[-37.3971023912332,-35.055603867752396,-20.15557608256847],[-35.90530048470914,-34.36180194986286,-20.480106343650313],[-34.419891792568286,-33.668057881359935,-20.777198203103005],[-32.94143261077337,-32.973973371837396,-21.045437860023057],[-31.47047923528728,-32.27915013088896,-21.28341151350693],[-30.00758796207291,-31.583189868108548,-21.489705362651087],[-28.553315087093033,-30.885694293089863,-21.66290560655199],[-27.10821690631051,-30.18626511542666,-21.80159844430617],[-25.672849715687498,-29.484504044713596,-21.90437007501248],[-24.247769811172624,-28.780012790561873,-21.96980669781984],[-22.833533488725926,-28.0723930625701,-21.996494511837835],[-21.430697044313952,-27.361246570328884,-21.9830197161521],[-20.039816773903247,-26.6461750234289,-21.9279685098483],[-18.661448973460445,-25.92678013146077,-21.829927092012078],[-17.29614993895222,-25.202663604015214,-21.68748166172912],[-15.94447596634516,-24.473427150682888,-21.49921841808508],[-14.606983351605841,-23.738672481054454,-21.2637235601656],[-13.284228390700711,-22.998001304720685,-20.97958328705691],[-11.976767379579904,-22.251015331284968,-20.645383797893654],[-10.685156614196224,-21.497316270350684,-20.259711291805676],[-9.409952390511855,-20.736505831514854,-19.821151967896615],[-8.151711004488838,-19.96818572437447,-19.328292025270006],[-6.910988752089359,-19.191957658526565,-18.77971766302947],[-5.688341929275552,-18.40742334356816,-18.174015080278586],[-4.484326832009543,-17.61418448909625,-17.509770476120988],[-3.299499756253475,-16.81184280470786,-16.78557004966025],[-2.1344169979695096,-15.999999999999995,-16.0],[21.1986654364892,0.2587549694507647,-0.2674088201578684]],"maxCircleResidualMM":0.0030623713756027736,"junctionContinuationMM":32.500681135661324,"id":1,"lengthMM":143.18731783393014,"turnRadians":1.6239255159898378,"nativeCenterChecks":16,"centersInsideSolid":0},{"face":102,"surface":"BSPLINE","radius":7.00005472565284,"points":[[-95.76202633310385,-65.89376104429411,-4.98168173379554e-12],[-93.32674258562909,-65.39805234922528,1.4599234622885224],[-92.00421384085806,-64.4069189642011,1.779486517270902],[-90.66693906389116,-63.430978463989504,2.1253356049150542],[-89.31547455069185,-62.46983255817987,2.496056924308522],[-87.95037659722374,-61.523082956360966,2.890236674538365],[-86.57220149945039,-60.59033136812155,3.3064610546915825],[-85.18150555333548,-59.67117950305032,3.743316263855266],[-83.77884505484252,-58.76522907073607,4.199388501116436],[-82.36477629993519,-57.87208178076761,4.673263965562151],[-80.93985558457712,-56.99133934273373,5.163528856279451],[-79.50463920472555,-56.12260346623524,5.668769372383238],[-78.05968345632931,-55.26547586089308,6.187571713031764],[-76.60554463535382,-54.41955823629959,6.718522077316493],[-75.14277903776531,-53.58445230204562,7.260206664325418],[-73.67194295952999,-52.75975976772201,7.811211673146575],[-72.19359269661402,-51.945082342919584,8.37012330286796],[-70.70828454498363,-51.14002173722921,8.935527752577595],[-69.21657480060503,-50.344179660241714,9.506011221363444],[-67.71901975944442,-49.557157821547975,10.080159908313567],[-66.21617571746185,-48.7785579307451,10.656560012535095],[-64.70859897059924,-48.0079816974487,11.23379773319061],[-63.196845814819106,-47.24503083125177,11.810459269375748],[-61.68147254608669,-46.48930704174498,12.385130820178682],[-60.16303546036715,-45.74041203851911,12.95639858468754],[-58.64209085362544,-44.997947531164655,13.522848761990542],[-57.119195021826734,-44.26151522927237,14.083067551175839],[-55.5949042609361,-43.5307168424328,14.635641151331596],[-54.06977486691873,-42.805154080236804,15.179155761545957],[-52.54436313573725,-42.08442865227979,15.712197580918806],[-51.01922536334238,-41.368142268183234,16.233352808610796],[-49.49491784569721,-40.65589663754282,16.74120764372258],[-47.97199687876708,-39.94729346994882,17.234348285342524],[-46.451018758517556,-39.2419344749916,17.711360932558947],[-44.93253978091419,-38.53942136226157,18.17083178446016],[-43.41711624192235,-37.83935584134904,18.61134704013444],[-41.90530443750763,-37.14133962184441,19.031492898670134],[-40.39766066363551,-36.44497441333801,19.42985555915556],[-38.894741216268805,-35.749861925422664,19.805021220685994],[-37.397102391345136,-35.05560386771526,20.15557608242451],[-35.90530048482399,-34.36180194981396,20.48010634347861],[-34.41989179267255,-33.66805788130968,20.77719820293535],[-32.94143261085815,-32.97397337179337,21.045437859881783],[-31.470479235347973,-32.27915013085597,21.283411513404964],[-30.007587962109277,-31.583189868088432,21.48970536259193],[-28.553315087109304,-30.885694293081738,21.662905606529755],[-27.108216906315356,-30.186265115426835,21.801598444305476],[-25.67284971569331,-29.48450404471574,21.904370075009744],[-24.247769811181392,-28.780012790564378,21.969806697809798],[-22.833533488734865,-28.072393062573436,21.996494511824228],[-21.430697044321242,-27.3612465703334,21.9830197161379],[-20.039816773907837,-26.64617502343455,21.92796850983567],[-18.66144897346208,-25.926780131467265,21.829927092002425],[-17.296149938951395,-25.202663604021932,21.687481661723034],[-15.944475966343155,-24.473427150688902,21.4992184180824],[-14.60698335160479,-23.73867248105854,21.26372356016536],[-13.284228390703426,-22.998001304721427,20.979583287057384],[-11.976767379585688,-22.25101533128389,20.645383797893768],[-10.685156614202363,-21.497316270350478,20.259711291805175],[-9.409952390516459,-20.73650583151732,19.821151967895357],[-8.151711004490842,-19.968185724380394,19.328292025268052],[-6.91098875208845,-19.191957658535713,18.779717663026997],[-5.68834192927217,-18.407423343579282,18.17401508027591],[-4.484326832004982,-17.6141844891072,17.509770476118558],[-3.299499756249787,-16.81184280471544,16.785570049658684],[-2.134416997969531,-15.999999999999993,15.999999999999996],[21.198665436322663,0.2587549695383018,0.2674088202513598]],"maxCircleResidualMM":0.0008925115149274365,"junctionContinuationMM":32.50068113554031,"id":2,"lengthMM":143.1873178337991,"turnRadians":1.6239254949308322,"nativeCenterChecks":16,"centersInsideSolid":0},{"face":125,"surface":"BSPLINE","radius":7.000054342734842,"points":[[-95.76202633310385,-65.89376104429411,-4.98168173379554e-12],[-98.19731008057502,-66.38946973935836,-1.718307886854421],[-97.31773537237069,-64.93034520935831,-2.1389982829573775],[-96.4152482479902,-63.398334217016526,-2.5982253214624293],[-95.49000282893803,-61.79722741138268,-3.093820886945725],[-94.54215323674455,-60.13081544168166,-3.6236168638897097],[-93.57185359294033,-58.40288895713862,-4.185445136776583],[-92.57925801905591,-56.61723860697861,-4.777137590088615],[-91.56452063662181,-54.77765504042663,-5.396526108308141],[-90.52779556716851,-52.887928906707806,-6.041442575917358],[-89.46923693222661,-50.95185085504717,-6.709718877398583],[-88.38899885332651,-48.973211534669424,-7.39918689723434],[-87.28723545198383,-46.95580159471788,-8.107678519961215],[-86.16410084971483,-44.903411684348896,-8.833025630105944],[-85.01974916804585,-42.8198324527779,-9.573060112155403],[-83.85433452850336,-40.7088545492204,-10.325613850596417],[-82.6680110526137,-38.57426862289177,-11.088518729915918],[-81.46093286190322,-36.41986532300747,-11.85960663460084],[-80.23325407789827,-34.24943529878294,-12.636709449138026],[-78.98512882212529,-32.06676919943357,-13.4176590580144],[-77.71671121611067,-29.87565767417436,-14.200287345717113],[-76.42815538139277,-27.67989137210781,-14.98242619679867],[-75.11961543951224,-25.48326094247199,-15.761907495739763],[-73.79124551196537,-23.289557034504803,-16.536563127009565],[-72.44319972028995,-21.102570297409777,-17.304224975118185],[-71.07563218600714,-18.926091380402,-18.06272492453166],[-69.68869703065445,-16.763910932620405,-18.80989485979678],[-68.28254837574372,-14.619819603270821,-19.543566665404374],[-66.85734034279099,-12.497608041605957,-20.261572225821162],[-65.41322705331248,-10.401066896878456,-20.961743425513895],[-63.950362628851785,-8.333986818240655,-21.641912149007286],[-62.46890119095392,-6.3001584548372085,-22.299910280830282],[-60.9689968611437,-4.303372455884982,-22.93356970546999],[-59.45080376094605,-2.3474194706010114,-23.540722307413404],[-57.91447601188589,-0.4360901482022333,-24.119199971147573],[-56.36016773548807,1.4268248620944084,-24.666834581159534],[-54.78803305327748,3.2375349110719447,-25.181458021936347],[-53.19822608677901,4.992249349513424,-25.660902177965063],[-51.59090095751754,6.687177528202019,-26.102998933732803],[-49.96621178699714,8.318528798029305,-26.505580173789156],[-48.32431269674151,9.882512509830955,-26.86647778265403],[-46.66535780829971,11.375338014330133,-27.183523644783605],[-44.98950124322109,12.793214662249884,-27.454549644633975],[-43.29689712303858,14.132351804365308,-27.67738766668731],[-41.58769956921474,15.38895879166991,-27.849869595535786],[-39.86206270328535,16.559244974913458,-27.96982731565059],[-38.12014064679895,17.639419704804514,-28.035092711482395],[-36.362087521303984,18.62569233205183,-28.043497667481954],[-34.58805744830993,19.514272207582284,-27.99287406821525],[-32.79820454932236,20.30136868229497,-27.881053798231072],[-30.992682945879103,20.983191106886885,-27.70586874197048],[-29.171646759518104,21.555948832055094,-27.46515078387458],[-27.335250111771607,22.015851208558804,-27.156731808418886],[-25.48364712414341,22.359107587441933,-26.778443700235137],[-23.61699191816073,22.581927319454852,-26.328118343790607],[-21.73543861535544,22.68051975529429,-25.803587623522628],[-19.83914133725944,22.651094245656974,-25.20268342386868],[-17.92825420539024,22.48986014145157,-24.523237629391662],[-16.0029313412711,22.193026793591777,-23.76308212465707],[-14.06332686644172,21.756803552796963,-22.920048794114955],[-12.109594902441764,21.177399769786444,-21.99196952221548],[-10.141889570803677,20.451024795309962,-20.976676193425657],[-8.160364993018682,19.57388798029045,-19.872000692310085],[-6.165175290603031,18.54219867554099,-18.67577490337365],[-4.156474585079167,17.35216623184852,-17.385830711106813],[-2.1344169979695025,15.999999999999995,-15.999999999999991],[21.494558400951952,0.19911340246692077,0.19427655698957302]],"maxCircleResidualMM":0.0006954273782175332,"junctionContinuationMM":32.714692247905944,"id":3,"lengthMM":184.42095652424905,"turnRadians":3.856241806873369,"nativeCenterChecks":16,"centersInsideSolid":0}],"bleedCount":0,"notes":"Internal swept/cylindrical centreline extraction. Merger and main-bore continuations are reduced tube unions, not a CAD-resolved flow field."},"nov7":{"source":"1'' Eductor 4 suction SHORT.step","sha256":"49263cce89b82f2d275148346730af0a94a94e372471c3a5de3c63b693002bc0","coordinates":"nozzle-relative CAD mm; flow +x","cadValid":true,"common":{"points":[[-114.70979579551081,-85.30330085889939,-2.4477358023231794e-16],[-95.76202633310388,-65.89376104429412,-4.974798351042864e-12]],"radius":7.0,"lengthMM":27.12467886289657},"branches":[{"face":98,"surface":"BSPLINE","radius":7.0000635705898056,"points":[[-95.76202633310388,-65.89376104429412,-4.974798351042864e-12],[-98.19731008059662,-66.38946973937996,1.718307886836651],[-97.31773537238976,-64.9303452093915,2.138998282949125],[-96.41524824800985,-63.398334217056956,2.598225321463055],[-95.49000282895807,-61.79722741142136,3.093820886947868],[-94.54215323676488,-60.13081544171296,3.623616863888014],[-93.57185359296085,-58.40288895715983,4.185445136768019],[-92.57925801907652,-56.61723860699016,4.777137590072327],[-91.56452063664236,-54.777655040432094,5.396526108285424],[-90.527795567189,-52.88792890671373,6.041442575891793],[-89.46923693224696,-50.951850855063256,6.7097188773759555],[-88.38899885334664,-48.973211534708405,7.3991868972225],[-87.28723545200589,-46.955801594777775,8.107678519959245],[-86.1641008497416,-44.90341168441732,8.833025630106569],[-85.0197491680787,-42.819832452845546,9.57306011215341],[-83.85433452854224,-40.708854549280844,10.32561385058872],[-82.66801105265705,-38.57426862294173,11.088518729901432],[-81.46093286194824,-36.41986532304655,11.85960663458054],[-80.23325407794059,-34.249435298813836,12.636709449114957],[-78.98512882215925,-32.066769199461895,13.417659057993674],[-77.71671121612908,-29.875657674208743,14.20028734570583],[-76.42815538139689,-27.679891372138133,14.982426196795723],[-75.11961543952552,-25.483260942500653,15.761907495731565],[-73.79124551199621,-23.28955703454072,16.53656312700099],[-72.44319972030158,-21.102570297437538,17.304224975105416],[-71.07563218603019,-18.92609138043513,18.062724924520325],[-69.6886970306821,-16.763910932660494,18.80989485978532],[-68.28254837577208,-14.619819603315191,19.54356666539389],[-66.85734034281603,-12.497608041646885,20.261572225811904],[-65.41322705332996,-10.401066896902957,20.961743425505404],[-63.9503626288654,-8.333986818246,21.641912149000436],[-62.4689011909717,-6.300158454831118,22.299910280826808],[-60.96899686117069,-4.303372455873824,22.93356970547075],[-59.45080376098415,-2.3474194705897866,23.540722307418378],[-57.914476011933864,-0.4360901481945942,24.119199971155876],[-56.36016773554166,1.4268248620962316,24.666834581169457],[-54.78803305332929,3.23753491106708,25.181458021945282],[-53.19822608681857,4.992249349502282,25.660902177969504],[-51.59090095753128,6.687177528186537,26.102998933728436],[-49.966211786991614,8.318528798013123,26.505580173779908],[-48.32431269673565,9.882512509817232,26.86647778264773],[-46.66535780830295,11.375338014319343,27.18352364478235],[-44.989501243232986,12.79321466224035,27.4545496446343],[-43.296897123052936,14.13235180435438,27.677387666683146],[-41.58769956923651,15.388958791660981,27.849869595532443],[-39.862062703314315,16.559244974908157,27.969827315651248],[-38.120140646826705,17.639419704801316,28.03509271148566],[-36.36208752131407,18.62569233204591,28.04349766748187],[-34.58805744831039,19.51427220757052,27.992874068211123],[-32.798204549335004,20.301368682279133,27.88105379822781],[-30.992682945905848,20.98319110687126,27.70586874197013],[-29.171646759540884,21.55594883204646,27.465150783876233],[-27.335250111763266,22.01585120856298,27.156731808418957],[-25.48364712412208,22.35910758744917,26.77844370023093],[-23.61699191814393,22.581927319458405,26.32811834378291],[-21.735438615351022,22.680519755293542,25.80358762351586],[-19.839141337265424,22.65109424565737,25.202683423870816],[-17.928254205399746,22.489860141456745,24.52323762940321],[-16.002931341279997,22.19302679359958,23.763082124670568],[-14.063326866448016,21.756803552804648,22.92004879412536],[-12.109594902445416,21.17739976979065,21.9919695222201],[-10.14188957080642,20.451024795305905,20.97667619342352],[-8.160364993021012,19.57388798027185,19.87200069229928],[-6.165175290605021,18.542198675511095,18.675774903356757],[-4.15647458508051,17.352166231821254,17.385830711091707],[-2.1344169979695238,15.999999999999991,16.000000000000004],[21.494558401330448,0.19911340254271792,-0.1942765570617233]],"maxCircleResidualMM":0.0007088845430032499,"junctionContinuationMM":32.71469224817845,"id":0,"lengthMM":184.42095652457198,"turnRadians":3.856241827940958,"nativeCenterChecks":16,"centersInsideSolid":0},{"face":99,"surface":"BSPLINE","radius":7.000282028264419,"points":[[-95.76202633310388,-65.89376104429412,-4.974798351042864e-12],[-93.32674258561471,-65.39805234921282,-1.4599234622906825],[-92.0042138408401,-64.40691896419946,-1.779486517280696],[-90.66693906386209,-63.43097846399638,-2.1253356049288743],[-89.31547455064768,-62.46983255819278,-2.4960569243234487],[-87.95037659716377,-61.523082956377316,-2.890236674552277],[-86.57220149937757,-60.59033136813875,-3.306461054703177],[-85.18150555325613,-59.67117950306581,-3.743316263863973],[-83.7788450547666,-58.765229070747296,-4.199388501122454],[-82.36477629987601,-57.87208178077186,-4.673263965566484],[-80.93985558455157,-56.99133934272834,-5.163528856283831],[-79.50463920474809,-56.12260346621677,-5.668769372382612],[-78.05968345640011,-55.265475860857215,-6.187571713006479],[-76.60554463546978,-54.4195582362445,-6.718522077254426],[-75.14277903792072,-53.58445230197206,-7.260206664223095],[-73.67194295971657,-52.75975976763336,-7.811211673009063],[-72.19359269682113,-51.94508234282193,-8.370123302708858],[-70.70828454519781,-51.14002173713114,-8.935527752419182],[-69.2165748008104,-50.3441796601545,-9.506011221236557],[-67.7190197596226,-49.5571578214855,-10.080159908257553],[-66.21617571759796,-48.778557930717554,-10.656560012578815],[-64.7085989707002,-48.00798169744414,-11.23379773329691],[-63.19684581489293,-47.24503083125876,-11.810459269508419],[-61.68147254613984,-46.489307041754856,-12.385130820309932],[-60.163035460404565,-45.740412038525875,-12.956398584798041],[-58.64209085365075,-44.99794753116528,-13.522848762069337],[-57.11919502184204,-44.26151522926653,-14.083067551220436],[-55.594904260942194,-43.530716842423146,-14.635641151347855],[-54.06977486691475,-42.80515408022855,-15.179155761548255],[-52.54436313572251,-42.08442865227853,-15.712197580923467],[-51.01922536332308,-41.36814226818104,-16.233352808602717],[-49.49491784567838,-40.655896637532315,-16.741207643688114],[-47.971996878751334,-39.94729346992616,-17.2343482852761],[-46.451018758504695,-39.24193447495634,-17.711360932463137],[-44.9325397809014,-38.539421362216665,-18.17083178434571],[-43.41711624190416,-37.83935584130087,-18.61134704002032],[-41.90530443747591,-37.14133962180277,-19.03149289858342],[-40.397660663579394,-36.444974413316125,-19.429855559131518],[-38.89474121617757,-35.74986192543475,-19.805021220761017],[-37.3971023912332,-35.055603867752396,-20.15557608256847],[-35.90530048470914,-34.36180194986286,-20.480106343650313],[-34.419891792568286,-33.668057881359935,-20.777198203103005],[-32.94143261077337,-32.973973371837396,-21.045437860023057],[-31.47047923528728,-32.27915013088896,-21.28341151350693],[-30.00758796207291,-31.583189868108548,-21.489705362651087],[-28.553315087093033,-30.885694293089863,-21.66290560655199],[-27.10821690631051,-30.18626511542666,-21.80159844430617],[-25.672849715687498,-29.484504044713596,-21.90437007501248],[-24.247769811172624,-28.780012790561873,-21.96980669781984],[-22.833533488725926,-28.0723930625701,-21.996494511837835],[-21.430697044313952,-27.361246570328884,-21.9830197161521],[-20.039816773903247,-26.6461750234289,-21.9279685098483],[-18.661448973460445,-25.92678013146077,-21.829927092012078],[-17.29614993895222,-25.202663604015214,-21.68748166172912],[-15.94447596634516,-24.473427150682888,-21.49921841808508],[-14.606983351605841,-23.738672481054454,-21.2637235601656],[-13.284228390700711,-22.998001304720685,-20.97958328705691],[-11.976767379579904,-22.251015331284968,-20.645383797893654],[-10.685156614196224,-21.497316270350684,-20.259711291805676],[-9.409952390511855,-20.736505831514854,-19.821151967896615],[-8.151711004488838,-19.96818572437447,-19.328292025270006],[-6.910988752089359,-19.191957658526565,-18.77971766302947],[-5.688341929275552,-18.40742334356816,-18.174015080278586],[-4.484326832009543,-17.61418448909625,-17.509770476120988],[-3.299499756253475,-16.81184280470786,-16.78557004966025],[-2.1344169979695096,-15.999999999999995,-16.0],[21.1986654364892,0.2587549694507647,-0.2674088201578684]],"maxCircleResidualMM":0.0030623713756027736,"junctionContinuationMM":32.500681135661324,"id":1,"lengthMM":143.1873178339302,"turnRadians":1.6239255159898338,"nativeCenterChecks":16,"centersInsideSolid":0},{"face":102,"surface":"BSPLINE","radius":7.00005472565284,"points":[[-95.76202633310388,-65.89376104429412,-4.974798351042864e-12],[-93.32674258562909,-65.39805234922528,1.4599234622885224],[-92.00421384085806,-64.4069189642011,1.779486517270902],[-90.66693906389116,-63.430978463989504,2.1253356049150542],[-89.31547455069185,-62.46983255817987,2.496056924308522],[-87.95037659722374,-61.523082956360966,2.890236674538365],[-86.57220149945039,-60.59033136812155,3.3064610546915825],[-85.18150555333548,-59.67117950305032,3.743316263855266],[-83.77884505484252,-58.76522907073607,4.199388501116436],[-82.36477629993519,-57.87208178076761,4.673263965562151],[-80.93985558457712,-56.99133934273373,5.163528856279451],[-79.50463920472555,-56.12260346623524,5.668769372383238],[-78.05968345632931,-55.26547586089308,6.187571713031764],[-76.60554463535382,-54.41955823629959,6.718522077316493],[-75.14277903776531,-53.58445230204562,7.260206664325418],[-73.67194295952999,-52.75975976772201,7.811211673146575],[-72.19359269661402,-51.945082342919584,8.37012330286796],[-70.70828454498363,-51.14002173722921,8.935527752577595],[-69.21657480060503,-50.344179660241714,9.506011221363444],[-67.71901975944442,-49.557157821547975,10.080159908313567],[-66.21617571746185,-48.7785579307451,10.656560012535095],[-64.70859897059924,-48.0079816974487,11.23379773319061],[-63.196845814819106,-47.24503083125177,11.810459269375748],[-61.68147254608669,-46.48930704174498,12.385130820178682],[-60.16303546036715,-45.74041203851911,12.95639858468754],[-58.64209085362544,-44.997947531164655,13.522848761990542],[-57.119195021826734,-44.26151522927237,14.083067551175839],[-55.5949042609361,-43.5307168424328,14.635641151331596],[-54.06977486691873,-42.805154080236804,15.179155761545957],[-52.54436313573725,-42.08442865227979,15.712197580918806],[-51.01922536334238,-41.368142268183234,16.233352808610796],[-49.49491784569721,-40.65589663754282,16.74120764372258],[-47.97199687876708,-39.94729346994882,17.234348285342524],[-46.451018758517556,-39.2419344749916,17.711360932558947],[-44.93253978091419,-38.53942136226157,18.17083178446016],[-43.41711624192235,-37.83935584134904,18.61134704013444],[-41.90530443750763,-37.14133962184441,19.031492898670134],[-40.39766066363551,-36.44497441333801,19.42985555915556],[-38.894741216268805,-35.749861925422664,19.805021220685994],[-37.397102391345136,-35.05560386771526,20.15557608242451],[-35.90530048482399,-34.36180194981396,20.48010634347861],[-34.41989179267255,-33.66805788130968,20.77719820293535],[-32.94143261085815,-32.97397337179337,21.045437859881783],[-31.470479235347973,-32.27915013085597,21.283411513404964],[-30.007587962109277,-31.583189868088432,21.48970536259193],[-28.553315087109304,-30.885694293081738,21.662905606529755],[-27.108216906315356,-30.186265115426835,21.801598444305476],[-25.67284971569331,-29.48450404471574,21.904370075009744],[-24.247769811181392,-28.780012790564378,21.969806697809798],[-22.833533488734865,-28.072393062573436,21.996494511824228],[-21.430697044321242,-27.3612465703334,21.9830197161379],[-20.039816773907837,-26.64617502343455,21.92796850983567],[-18.66144897346208,-25.926780131467265,21.829927092002425],[-17.296149938951395,-25.202663604021932,21.687481661723034],[-15.944475966343155,-24.473427150688902,21.4992184180824],[-14.60698335160479,-23.73867248105854,21.26372356016536],[-13.284228390703426,-22.998001304721427,20.979583287057384],[-11.976767379585688,-22.25101533128389,20.645383797893768],[-10.685156614202363,-21.497316270350478,20.259711291805175],[-9.409952390516459,-20.73650583151732,19.821151967895357],[-8.151711004490842,-19.968185724380394,19.328292025268052],[-6.91098875208845,-19.191957658535713,18.779717663026997],[-5.68834192927217,-18.407423343579282,18.17401508027591],[-4.484326832004982,-17.6141844891072,17.509770476118558],[-3.299499756249787,-16.81184280471544,16.785570049658684],[-2.134416997969531,-15.999999999999993,15.999999999999996],[21.198665436322663,0.2587549695383018,0.2674088202513598]],"maxCircleResidualMM":0.0008925115149274365,"junctionContinuationMM":32.50068113554031,"id":2,"lengthMM":143.1873178337991,"turnRadians":1.6239254949308255,"nativeCenterChecks":16,"centersInsideSolid":0},{"face":118,"surface":"BSPLINE","radius":7.000054342734842,"points":[[-95.76202633310388,-65.89376104429412,-4.974798351042864e-12],[-98.19731008057508,-66.38946973935843,-1.7183078868543902],[-97.31773537237063,-64.93034520935835,-2.138998282957378],[-96.4152482479902,-63.39833421701653,-2.598225321462419],[-95.490002828938,-61.797227411382686,-3.093820886945729],[-94.54215323674453,-60.130815441681676,-3.623616863889716],[-93.5718535929403,-58.402888957138615,-4.185445136776587],[-92.57925801905591,-56.61723860697862,-4.777137590088611],[-91.56452063662178,-54.777655040426616,-5.396526108308141],[-90.52779556716851,-52.8879289067078,-6.041442575917369],[-89.46923693222661,-50.95185085504718,-6.709718877398585],[-88.38899885332651,-48.973211534669424,-7.399186897234339],[-87.28723545198383,-46.95580159471788,-8.107678519961215],[-86.16410084971483,-44.90341168434888,-8.833025630105956],[-85.01974916804585,-42.81983245277791,-9.573060112155392],[-83.85433452850339,-40.7088545492204,-10.325613850596415],[-82.66801105261368,-38.574268622891786,-11.08851872991592],[-81.46093286190319,-36.41986532300748,-11.859606634600844],[-80.23325407789827,-34.249435298782934,-12.636709449138026],[-78.98512882212529,-32.06676919943357,-13.417659058014397],[-77.71671121611064,-29.87565767417435,-14.200287345717122],[-76.42815538139277,-27.67989137210781,-14.98242619679867],[-75.11961543951222,-25.483260942471997,-15.761907495739761],[-73.79124551196537,-23.289557034504803,-16.53656312700956],[-72.44319972028995,-21.102570297409784,-17.30422497511817],[-71.07563218600714,-18.926091380402,-18.06272492453166],[-69.68869703065445,-16.76391093262041,-18.809894859796785],[-68.28254837574372,-14.619819603270829,-19.543566665404374],[-66.85734034279099,-12.497608041605957,-20.261572225821165],[-65.41322705331248,-10.40106689687846,-20.961743425513895],[-63.950362628851785,-8.333986818240659,-21.64191214900729],[-62.46890119095392,-6.3001584548372165,-22.299910280830275],[-60.9689968611437,-4.303372455884986,-22.933569705469985],[-59.450803760946066,-2.347419470601005,-23.5407223074134],[-57.9144760118859,-0.4360901482022346,-24.119199971147562],[-56.36016773548807,1.4268248620944095,-24.666834581159534],[-54.78803305327748,3.23753491107194,-25.181458021936347],[-53.19822608677901,4.992249349513424,-25.660902177965063],[-51.59090095751754,6.687177528202029,-26.102998933732795],[-49.96621178699714,8.318528798029309,-26.50558017378916],[-48.324312696741494,9.882512509830942,-26.86647778265403],[-46.665357808299696,11.375338014330142,-27.183523644783605],[-44.98950124322109,12.793214662249882,-27.454549644633975],[-43.2968971230386,14.132351804365321,-27.677387666687302],[-41.58769956921477,15.388958791669909,-27.849869595535786],[-39.86206270328535,16.55924497491346,-27.96982731565059],[-38.120140646798966,17.639419704804514,-28.035092711482395],[-36.362087521303984,18.625692332051834,-28.043497667481958],[-34.58805744830993,19.51427220758229,-27.992874068215254],[-32.79820454932236,20.30136868229497,-27.88105379823107],[-30.99268294587909,20.983191106886903,-27.705868741970473],[-29.171646759518104,21.555948832055094,-27.46515078387458],[-27.335250111771607,22.01585120855881,-27.156731808418886],[-25.48364712414341,22.359107587441926,-26.778443700235133],[-23.616991918160743,22.581927319454852,-26.32811834379061],[-21.73543861535544,22.680519755294295,-25.803587623522624],[-19.83914133725944,22.651094245656992,-25.202683423868677],[-17.92825420539024,22.489860141451555,-24.523237629391666],[-16.0029313412711,22.19302679359178,-23.763082124657064],[-14.06332686644172,21.75680355279696,-22.92004879411496],[-12.109594902441778,21.177399769786447,-21.991969522215495],[-10.141889570803677,20.451024795309962,-20.976676193425657],[-8.160364993018682,19.57388798029045,-19.872000692310085],[-6.165175290603031,18.54219867554099,-18.675774903373654],[-4.15647458507916,17.352166231848525,-17.385830711106806],[-2.1344169979695025,15.99999999999999,-15.999999999999986],[21.494558400951796,0.19911340246686215,0.1942765569895073]],"maxCircleResidualMM":0.0006954273782163134,"junctionContinuationMM":32.71469224790582,"id":3,"lengthMM":184.42095652424905,"turnRadians":3.856241806873207,"nativeCenterChecks":16,"centersInsideSolid":0}],"bleedCount":0,"notes":"Internal swept/cylindrical centreline extraction. Merger and main-bore continuations are reduced tube unions, not a CAD-resolved flow field."}}
;
const EDUCTOR_CALIBRATION={
  "version": 2,
  "scope": "Historical British body only (revC); 3.5 mm calibration, 2.5 mm holdout. Not transferred to 20Nov CAD.",
  "trainingCase": "british-p3bar-sst-v7",
  "trainingHash": "cfa61df0b09ba35b0b9be170abb80332416b2869def9b1cf6375ce97fea45a83",
  "dataHash": "fa8b102b76bcc0d85cecaaae429f33e4e26e3fee00f33a91bfbaf5bf48abb012",
  "coefficients": {
    "feederResistance": 4.197021050729741,
    "mixingLoss": 0.7413255389701133,
    "turbulenceScale": 0.72334973667736,
    "spreadScale": 0.9497442890891235,
    "developmentD": 1.3397478822494273
  },
  "targets": {
    "flows": {
      "motive_g_s": 6.643176784387476,
      "feeder_g_s": 5.7894359934826,
      "bleed_g_s": 7.680958758302698,
      "outlet_g_s": 20.116178967675932,
      "inlet_sum_g_s": 20.113571536172774,
      "imbalance_pct": 0.012963543040919533
    },
    "peak": {
      "G_s_inv": 3829961.607121514,
      "position_m": [
        0.0071295211617202125,
        0.0018790282577896128,
        -0.0006263191045274916
      ],
      "x_D": 2.037006046205775
    },
    "shearRmsRadius": 0.00265459375756306
  },
  "flowFit": {
    "best": [
      1.4343750000000002,
      -0.2993154263142097
    ],
    "score": 1.1172669398056097e-7,
    "trace": [
      {
        "round": 0,
        "score": 0.0034573262527046586,
        "parameters": [
          1.2,
          -0.2231435513142097
        ]
      },
      {
        "round": 1,
        "score": 0.001174923232386583,
        "parameters": [
          1.5,
          -0.2231435513142097
        ]
      },
      {
        "round": 2,
        "score": 0.0004519762484833406,
        "parameters": [
          1.5,
          -0.3481435513142097
        ]
      },
      {
        "round": 3,
        "score": 0.000020064609015007302,
        "parameters": [
          1.425,
          -0.2856435513142097
        ]
      },
      {
        "round": 4,
        "score": 0.000020064609015007302,
        "parameters": [
          1.425,
          -0.2856435513142097
        ]
      },
      {
        "round": 5,
        "score": 0.000006412256934044864,
        "parameters": [
          1.44375,
          -0.3012685513142097
        ]
      },
      {
        "round": 6,
        "score": 9.639152391928985e-7,
        "parameters": [
          1.4343750000000002,
          -0.3012685513142097
        ]
      },
      {
        "round": 7,
        "score": 1.4488517963372385e-7,
        "parameters": [
          1.4343750000000002,
          -0.2973623013142097
        ]
      },
      {
        "round": 8,
        "score": 1.1172669398056097e-7,
        "parameters": [
          1.4343750000000002,
          -0.2993154263142097
        ]
      }
    ]
  },
  "fieldFit": {
    "best": [
      -0.32386244393873237,
      -0.051562500000000046,
      0.29248144868579035
    ],
    "score": 0.0008133029425385216,
    "trace": [
      {
        "round": 0,
        "score": 0.010573628078886423,
        "parameters": [
          -0.35667494393873245,
          -0.3,
          0.37685644868579027
        ]
      },
      {
        "round": 1,
        "score": 0.003755152818286751,
        "parameters": [
          -0.35667494393873245,
          -0.15,
          0.37685644868579027
        ]
      },
      {
        "round": 2,
        "score": 0.0019397353379613362,
        "parameters": [
          -0.35667494393873245,
          -0.15,
          0.30185644868579026
        ]
      },
      {
        "round": 3,
        "score": 0.001808047517079651,
        "parameters": [
          -0.35667494393873245,
          -0.15,
          0.2643564486857903
        ]
      },
      {
        "round": 4,
        "score": 0.0015586466371131848,
        "parameters": [
          -0.35667494393873245,
          -0.13125,
          0.2643564486857903
        ]
      },
      {
        "round": 5,
        "score": 0.0008775483432170918,
        "parameters": [
          -0.3285499439387324,
          -0.05625000000000004,
          0.29248144868579035
        ]
      },
      {
        "round": 6,
        "score": 0.0008133029425385216,
        "parameters": [
          -0.32386244393873237,
          -0.051562500000000046,
          0.29248144868579035
        ]
      }
    ]
  },
  "limitations": "A fit is not independent validation. Losses are effective coefficients, not measured component K values. The finite mixing-length closure and local scalar distribution are not CFD-resolved."
}
;
'use strict';
/* Revision 2: explicit suction circuit, separated bleed/mixing control volumes.
 * SI internally. Closure coefficients are assumptions, never CFD observations.
 */
const EductorUpgrades=(()=>{
 const CONT=typeof EductorContinuation!=='undefined'?EductorContinuation:require('./continuation.js');
 const N=typeof EductorNetwork!=='undefined'?EductorNetwork:require('../exposure/network.js');
 const CAD=typeof EDUCTOR_TOPOLOGY!=='undefined'?EDUCTOR_TOPOLOGY:require('./topology.json');
 const candidate=typeof EDUCTOR_CANDIDATE!=='undefined'?EDUCTOR_CANDIDATE:require('../phases/calibration-candidate.json');
 const calibration=typeof EDUCTOR_CALIBRATION!=='undefined'?EDUCTOR_CALIBRATION:require('./calibration.json');
 function calibratedGeometry(g,p){if(p.calibrationCandidate){const c=candidate.fits[g.family];if(c)return {...g,feederResistance:c.feederResistance,mixingLoss:c.mixingLoss,calibrationCase:'Phase 3 family candidate; held-out validation failed'};}return g.key==='revC'&&p.calibration!==false?{...g,...calibration.coefficients,calibrationCase:calibration.trainingCase}:g;}
 const area=d=>Math.PI*d*d/4,clamp=(v,a,b)=>Math.max(a,Math.min(b,v)),norm=v=>Math.hypot(...v);
 const sub=(a,b)=>a.map((v,i)=>v-b[i]),dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
 function topology(g,release=30){
  const source=CAD[g.key||('nov'+g.id)];
  if(source)return source;
  const sn=Math.sin(g.angle*Math.PI/180),cs=Math.cos(g.angle*Math.PI/180),cap=Math.max(release+8,g.db*1.5);
  return {coordinates:'nozzle-relative CAD mm; flow +x',cadValid:null,source:'Nominal geometry, no matching STEP',common:{points:[[-cs*cap,0,sn*cap],[-cs*cap,0,sn*cap]],radius:(g.df||g.db)/2,lengthMM:0},branches:[{id:0,radius:g.db/2,points:[[-cs*cap,0,sn*cap],[0,0,0]],lengthMM:cap,turnRadians:g.angle*Math.PI/180}],notes:'Nominal straight suction passage; CAD topology unavailable for this historical geometry.'};
 }
 function friction(Re,roughness=0){
  if(!(Re>0))return 0;
  const lam=64/Re,turb=(-1.8*Math.log10((roughness/3.7)**1.11+6.9/Math.max(Re,1)))**-2;
  return Re<2300?lam:Re>4000?turb:lam+(turb-lam)*(Re-2300)/1700;
 }
 function hydraulic(g,P,Ps,T,md,loss=1){
  const topo=topology(g),visc=1.81e-5,dp=Math.max(0,Ps-P),rho=Ps/(N.R*T);
  let guess=Math.max(1e-8,area(g.db/1000)*Math.sqrt(2*rho*dp)),branches,common,Aeff,s;
  for(let iter=0;iter<8;iter++){
   common={...topo.common,A:area(topo.common.radius*2/1000)};
   const Re=guess*(2*common.radius/1000)/(visc*common.A);
   common.K=loss*(.5+friction(Re)*common.lengthMM/(2*common.radius))*(g.feederResistance??1);
   branches=topo.branches.map(b=>{const A=area(2*b.radius/1000),re=guess/topo.branches.length*(2*b.radius/1000)/(visc*A),K=loss*(.5+friction(re)*b.lengthMM/(2*b.radius)+.2*b.turnRadians);return {...b,A,K,conductance:A/Math.sqrt(Math.max(K,.05))};});
   const conductance=branches.reduce((s,b)=>s+b.conductance,0);
   Aeff=1/Math.sqrt(common.K/common.A**2+1/conductance**2);
   s=N.inlet(Ps,P,T,Aeff,md);guess=.5*(guess+s.m);
  }
  const C=branches.reduce((v,b)=>v+b.conductance,0),totalResistance=common.K/common.A**2+1/C**2;
  common.pressureDropFraction=(common.K/common.A**2)/totalResistance;
  const junctionP=Math.sqrt(Math.max(P*P,Ps*Ps-common.pressureDropFraction*(Ps*Ps-P*P)));
  branches=branches.map(b=>({...b,fraction:b.conductance/C,m:s.m*b.conductance/C,powder:md*b.conductance/C,Re:s.m*b.conductance/C*(2*b.radius/1000)/(visc*b.A),pressureIn:junctionP,pressureOut:P}));
  return {...s,topology:topo,branches,common:{...common,pressureIn:Ps,pressureOut:junctionP},effectiveArea:Aeff,closure:'Darcy–Weisbach series/parallel resistance; isentropic restriction for choking. Junction and bend losses remain uncalibrated.'};
 }
 function holes(g,P,Po,T,loss){
  const count=g.family==='british'?(g.nh||8):0,each=area((g.dh||8)/1000)/Math.sqrt(1+1.5*loss),ports=[];
  for(let i=0;i<count;i++)ports.push({...N.inlet(Po,P,T,each),id:i,angle:2*Math.PI*i/count});
  const h=ports[0]||{m:0,u:0,T,T0:T,p:P,p0:Po,A:0,M:0,choked:false};
  return {...h,m:h.m*count,A:each*count,ports};
 }
 function train(g,p){
  const Po=p.ps*1e5,Tm=p.temp+273.15,md=p.feed/3600,loss=p.loss??1;
  const geometries=[{key:'revB1',dn:3.5,db:34.8,d:34.8,de:34.8,L:110.3,angle:90,nports:1,family:'branched'},{key:'revB2',dn:3.5,db:22.1,d:22.1,de:22.1,L:61.7+(p.ext||0)*1000,angle:45,nports:1,family:'branched'}];
  const zero={m:0,u:0,T0:293.15,T:293.15,p:Po,p0:Po,A:0,M:0,choked:false};
  function stage(gg,P,p0,Ps,Ts,incoming=null){
   const d=N.dimensions(gg,loss),n=p0===null?{...zero,p0:Ps,T0:Tm}:N.inlet(p0,P,Tm,d.An),s=incoming||hydraulic(gg,P,Ps,Ts,md,loss);
   if(!(s.m>0))return null;
   const mix=N.mixing(P,d.A,n,s,zero,md,gg.angle,d.K+(gg.mixingLoss??.8)*loss);if(!mix)return null;
   const out=N.duct(mix.m,md,mix.T0,mix.p0,d.Ae);return out?{g:gg,P,n,s,h:zero,mix,out,d}:null;
  }
  const result=N.solvePressureRoot(P=>{
   const a=stage(geometries[0],P,p.off1?null:(p.p01??p.p0)*1e5,Po-(p.feedDrop||0),293.15);if(!a)return null;
   const pt=a.out.p0-(p.interDrop||0),Ts=a.mix.T0;
   const incoming=N.solvePressureRoot(P2=>{const s=hydraulic(geometries[1],P2,pt,Ts,md,loss);return {P:P2,s,error:(s.m-a.mix.m)*1e7};},pt);if(!incoming)return null;
   const b=stage(geometries[1],incoming.P,p.p0*1e5,pt,Ts,incoming.s);if(!b)return null;
   return {stages:[a,b],error:b.out.p-Po};},Math.min(Po-(p.feedDrop||0),(p.p01??p.p0)*1e5));
  if(!result)return {ok:false,attempted:true,constraints:['coupled-staged-outlet-mismatch'],reason:'The shared staged train closure has no admissible coupled outlet root at these settings.'};
  const stages=result.stages.map(a=>{
   const {n,s,h,mix,out,P}=a,thrust=n.m*n.u+(n.p-P)*n.A;
   return {ok:true,modelVersion:2,manual:false,caution:true,compressible:Math.max(mix.M,out.M,s.M)>.3,closure:'coupled-staged-network-v3',reason:'Both stages use the same hydraulic circuit and mixing-loss closure. Interstage gas mass, total temperature and total pressure are transferred.',mp:n.m,ms:s.m,mh:0,mg:mix.m,md,p0:n.p0,ps:Po,po:out.p,systemOutletP:Po,inletP:s.p0,Tin:s.T0,T0:n.T0,T:out.T,totalT:mix.T0,portP:P,throatP:mix.p,uJet:n.m?Math.sqrt(2*N.cp*n.T0*(1-(P/n.p0)**((N.gamma-1)/N.gamma))):0,power:n.m*N.cp*n.T0*(1-(P/n.p0)**((N.gamma-1)/N.gamma))+s.m*N.cp*s.T0*(1-(P/s.p0)**((N.gamma-1)/N.gamma)),choked:n.choked,circuit:{...s,flowUsed:true},gate:{needed:Math.max(mix.M,out.M,s.M)>.3,reason:'Shared staged train balance'},network:{regime:'Coupled staged train',bulkMach:Math.max(mix.M,out.M,s.M),mixingMach:mix.M,outletMach:out.M,pressureResidual:result.error,residuals:mix.residuals,entropyProduction:mix.entropy,motiveThrust:thrust,mixingPressure:mix.p,mixingTemperature:mix.T,mixingVelocity:mix.u,outletVelocity:out.u,secondaryVelocity:s.u,bleedVelocity:0,assumptions:'Same staged hydraulic and mixing-loss closure in both stages. Isentropic transfer minus specified interstage pressure loss. Particle releases are stage-local.'}};
  });
  return {...stages[g.key==='revB2'?1:0],stages};
 }
 function solve(g,p){
  g=calibratedGeometry(g,p);
  const Po=p.ps*1e5,Ps=Po-(p.feedDrop||p.interDrop||0),p0=p.p0*1e5,Tm=p.temp+273.15,md=p.feed/3600,loss=p.loss??1;
  if(!(p0>Po&&Po>0&&Ps>0&&Tm>0&&md>=0&&loss>0))return {ok:false,reason:'Motive pressure must exceed the shared positive vessel/outlet pressure.'};
  if(g.key?.startsWith('revB')&&p.trainClosure==='staged')return train(g,p);
  // Retain the previous train as an explicitly labelled comparison closure.
  if(g.key?.startsWith('revB')){const f=N.solve(g,p);if(f.ok){f.modelVersion=2;f.upgradeScope='Coupled conservative train; nominal stage-local passages';f.circuit=hydraulic(g,f.portP,f.inletP,f.Tin,md,loss);f.gate={needed:f.network.bulkMach>=.3,reason:'Coupled train compressible energy and momentum balance'};}return f;}
  const dimensions=N.dimensions(g,loss),Ah=dimensions.Ah,zero={m:0,u:0,T:293.15,T0:293.15,p:Po,p0:Po,A:0,M:0,choked:false};
  const constraints=new Set();
  function evaluate(P){
   const n=N.inlet(p0,P,Tm,dimensions.An),s=hydraulic(g,P,Ps,293.15,md,loss);
   if(!(n.m>0&&s.m>1e-10)){constraints.add('forward-inlet-flow');return null;}
   let h=zero,mix,out,pre=null,bleedPressure=P;
   // Pressure recovery between the suction section and the radial holes.
   // This is a finite mixing-length closure, not instantaneous complete mixing
   // upstream of the holes (which incorrectly predicts zero bleed entrainment).
   const bleedX=(g.xh!=null?Math.max(0,(g.xh-(g.xn||0))/1000):0),mixLength=8*g.dn/1000;
   bleedPressure=P+(Po-P)*(Ah?1-Math.exp(-bleedX/mixLength):0);
   h=holes(g,bleedPressure,Po,293.15,loss);
   const mixingLoss=(g.mixingLoss??.8)*loss;
   const diagnostic={};mix=N.mixing(P,dimensions.A,n,s,h,md,g.angle??45,dimensions.K+mixingLoss,diagnostic);
   if(!mix){constraints.add(diagnostic.constraint);return null;}
   out=N.duct(mix.m,md,mix.T0,mix.p0,dimensions.Ae);
   if(!out){constraints.add('outlet-choking');return null;}
   return {P,n,s,h,mix,out,pre,bleedPressure,error:out.p-Po};
  }
  if(p.returnEvaluator)return {evaluate,limit:Math.min(Ps,p0)};
  let a=N.solvePressureRoot(evaluate,Math.min(Ps,p0),{seeds:p.pressureSeed?[p.pressureSeed]:[]});
  let continuation=null;
  if(!a&&p.continuation!==false){
   // Deterministic neighbours, never a cached result from another UI run.
   for(const direction of [-1,1]){
    const near=[.98,.96].map(scale=>solve(g,{...p,ps:p.ps*(direction<0?scale:2-scale),continuation:false}));
    if(!near.every(q=>q.ok&&!q.network.shock))continue;
    const a0=[near[1].portP/1e5,p.ps*(direction<0?.96:1.04)],b0=[near[0].portP/1e5,p.ps*(direction<0?.98:1.02)];
    const residual=(x,y)=>{if(!(x>0&&y>0&&y<p.p0))return NaN;const q=solve(g,{...p,ps:y,returnEvaluator:true,continuation:false}).evaluate(x*1e5);return q?q.error/1e5:NaN;};
    const path=CONT.trace(residual,a0,b0,{step:.01*p.ps,count:16});
    const seeds=path.filter((q,i)=>i&&((path[i-1][1]-p.ps)*(q[1]-p.ps)<=0)).map((q,i)=>q[0]*1e5);
    a=N.solvePressureRoot(evaluate,Math.min(Ps,p0),{seeds,dense:true});continuation={method:'pseudo-arclength',points:path};if(a)break;
   }
  }
  // Solve the sonic mixing / internal-shock branch using the SAME inlet
  // circuit and bleed recovery closure, without reverting to lumped flows.
  if(!a){
   if(dimensions.Ae>dimensions.A){
    const critical=N.solvePressureRoot(P=>{
     const n=N.inlet(p0,P,Tm,dimensions.An),s=hydraulic(g,P,Ps,293.15,md,loss);
     if(!n.choked||s.m<=1e-10)return null;
     const recovery=1-Math.exp(-(g.xh!=null?Math.max(0,g.xh-(g.xn||0)):0)/(8*g.dn)),bleedPressure=P+(Po-P)*recovery,h=holes(g,bleedPressure,Po,293.15,loss);
     const m=n.m+s.m+h.m,mt=m+md,T0=(n.m*Tm+(s.m+h.m)*293.15)/m,T=T0/1.2,u=Math.sqrt(N.gamma*N.R*T/(1+md/m)),pressure=m*N.R*T/(dimensions.A*u),totalP=pressure*(T0/T)**3.5,wall=(dimensions.K+(g.mixingLoss??.8)*loss)/2*mt*u;
     const I=P*dimensions.A+n.m*n.u+(n.p-P)*n.A+((s.m+md)*s.u+(s.p-P)*s.A)*Math.cos((g.angle??90)*Math.PI/180),residual=pressure*dimensions.A+mt*u+wall-I;
     const entropy=m*(N.cp*Math.log(T0)-N.R*Math.log(totalP))-[n,s,h].reduce((v,q)=>v+q.m*(N.cp*Math.log(q.T0)-N.R*Math.log(q.p0)),0);
     if(entropy<0)return null;
     const mix={m,md,T0,p0:totalP,A:dimensions.A,u,T,p:pressure,rho:pressure/(N.R*T),M:1,wall,entropy,residuals:{mass:m-pressure/(N.R*T)*dimensions.A*u,momentum:residual,energy:m*N.cp*T+.5*mt*u*u-m*N.cp*T0}};
     return {P,n,s,h,mix,bleedPressure,error:residual/dimensions.A};
    },Math.min(Ps,p0),{dense:true,accept:q=>!!N.shockDiffuser(q.mix,dimensions.Ae,Po)});
    const shock=critical&&N.shockDiffuser(critical.mix,dimensions.Ae,Po);
    if(shock)a={...critical,shock,out:shock.out,error:shock.error};
   }
   if(!a){constraints.add('critical-branch-or-outlet-mismatch');return {ok:false,attempted:true,modelVersion:2,constraints:[...constraints],continuation,reason:'No admissible branch. Constraints: '+[...constraints].join(', ')+'. Subsonic, continuation and critical diffuser branches checked; no substitute airflow.'};}
  }
  const {n,s,h,mix,out}=a;
  const power=n.m*N.cp*Tm*(1-(a.P/p0)**((N.gamma-1)/N.gamma))+s.m*N.cp*s.T0*(1-(a.P/Ps)**((N.gamma-1)/N.gamma))+h.m*N.cp*293.15*(1-(a.bleedPressure/Po)**((N.gamma-1)/N.gamma));
  const bulkMach=Math.max(mix.M,out.M,s.M,h.M,a.shock?.before.M||0),compressible=bulkMach>=.3||(Ps-a.P)/Ps>.05;
  const network={continuation,regime:compressible?'Staged compressible balance':'Staged low-Mach balance',staged:true,bulkMach,mixingMach:mix.M,outletMach:out.M,pressureResidual:a.error,residuals:mix.residuals,upstreamResiduals:a.pre?.residuals,entropyProduction:mix.entropy+(a.pre?.entropy||0),mixingPressure:mix.p,mixingTemperature:mix.T,mixingVelocity:mix.u,outletVelocity:out.u,secondaryVelocity:s.u,bleedVelocity:h.u,motiveThrust:n.m*n.u+(n.p-a.P)*n.A,motiveChoked:n.choked,secondaryChoked:s.choked,bleedChoked:h.choked,bleedPressure:a.bleedPressure,assumptions:'Separate suction and radial-hole pressures with finite mixing-length pressure recovery. Global mass, axial momentum and total energy balance includes a mixing/wall loss. Outlet area change is isentropic. This does not resolve shock cells or recirculation; coefficients require geometry-specific validation.'};
  if(a.shock){network.shock={...a.shock,diameter:Math.sqrt(4*a.shock.area/Math.PI)};network.regime='Critical mixing · internal diffuser shock';network.assumptions+=' One ideal normal shock matches the outlet; shock-particle interaction is unresolved.';}
  return {ok:true,modelVersion:2,calibration:g.calibrationCase||null,manual:false,compressible,caution:true,closure:'staged-network-v2',reason:'CAD/nominal suction circuit with staged pressure recovery and individual bleed restrictions.',gate:{needed:compressible,reason:compressible?'Bulk Mach or density change requires compressible evaluation':'Bulk Mach < 0.3 and density variation < 5%; low-Mach regime'},mp:n.m,ms:s.m,mh:h.m,mg:mix.m,md,p0,ps:Po,po:out.p,systemOutletP:Po,inletP:Ps,Tin:293.15,T0:Tm,T:out.T,totalT:mix.T0,portP:a.P,throatP:mix.p,uJet:Math.sqrt(2*N.cp*Tm*(1-(a.P/p0)**((N.gamma-1)/N.gamma))),power,choked:n.choked,network,circuit:{...s,flowUsed:true},bleedPorts:h.ports||[]};
 }
 return {topology,friction,hydraulic,holes,solve,calibratedGeometry,calibration,CAD,area,clamp,norm,sub,dot};
})();
if(typeof module!=='undefined')module.exports=EductorUpgrades;

'use strict';
/* CAD-centreline 3D tube union + conservative dilute scalar reconstruction.
 * The internal junctions and turbulent field remain reduced approximations.
 */
const EductorSpatial=(()=>{
 const WALL=typeof EductorWalls!=='undefined'?EductorWalls:require('./walls.js');
 const JET=typeof EductorJet!=='undefined'?EductorJet:require('./jet.js');
 const NUM=typeof ExposureNumerics!=='undefined'?ExposureNumerics:require('../exposure/numerics.js');
 const U=typeof EductorUpgrades!=='undefined'?EductorUpgrades:require('./physics.js');
 const {clamp,area}=U,R=287.05,cp=1004.675;
 const smooth=x=>(x=clamp(x,0,1),x*x*(3-2*x));
 const mu=T=>1.716e-5*(T/273.15)**1.5*(273.15+110.4)/(T+110.4);
 function geometry(g,p,legacy){
  const source=U.topology(g,p.release),cad=!!U.CAD[g.key],convert=q=>[q[0]/1000,q[2]/1000,q[1]/1000];
  const branches=source.branches.map((b,i)=>{
   let pts=b.points.filter((v,j)=>j===0||j===b.points.length-1||j%3===0).map(convert);
   if(source.common.lengthMM>0)pts=[convert(source.common.points[0]),...pts];
   const segments=[];let length=0;
   for(let j=0;j<pts.length-1;j++){const a=pts[j],b=pts[j+1],delta=b.map((v,k)=>v-a[k]),L=Math.hypot(...delta);if(L<1e-10)continue;segments.push({a,b,L,s:length,tangent:delta.map(v=>v/L),radius:(j===0&&source.common.lengthMM>0?source.common.radius:source.branches[i].radius)/1000});length+=L;}
   const last=segments.at(-1),entry=last.b;
   return {id:i,segments,length,entry,angle:Math.atan2(last.a[2],last.a[1]),radius:b.radius/1000};
  });
  function suctionRadius(y){const pts=g.suctionProfile;if(!pts||p.flowPhysics!=='phase2')return null;const mm=clamp(y*1000,pts[0][0],pts.at(-1)[0]);let i=1;while(i<pts.length-1&&pts[i][0]<mm)i++;const a=pts[i-1],b=pts[i],slope=(b[1]-a[1])/(b[0]-a[0]);return {radius:(a[1]+slope*(mm-a[0]))/1000,slope:y*1000<=pts[0][0]||y*1000>=pts.at(-1)[0]?0:slope};}
  function closest(x,y,z=0,only=null){
   let best=null,bestClear=-Infinity;
   for(let bi=0;bi<branches.length;bi++){if(only!==null&&bi!==only)continue;const branch=branches[bi];for(const segment of branch.segments){
    const a=segment.a,v=segment.tangent,t=clamp((x-a[0])*v[0]+(y-a[1])*v[1]+(z-a[2])*v[2],0,segment.L),cx=a[0]+t*v[0],cy=a[1]+t*v[1],cz=a[2]+t*v[2],distance=Math.hypot(x-cx,y-cy,z-cz),shape=suctionRadius(cy),radius=shape?.radius??segment.radius,clearance=radius-distance;
    if(clearance>bestClear){bestClear=clearance;best={branch,segment:shape?{...segment,radius}:segment,radiusSlope:shape?shape.slope*segment.tangent[1]:0,t,cx,cy,cz,distance,clearance,s:segment.s+t};}
   }}
   best.center=[best.cx,best.cy,best.cz];best.offset=[x-best.cx,y-best.cy,z-best.cz];return best;
  }
  function clearance(x,y,dp=0,z=0){
   const rr=legacy.radius(x),main=Math.min(x,legacy.end-x,rr-Math.hypot(y,z));
   if(main>dp/2+1e-8)return main-dp/2;
   let tube=-Infinity;
   for(const b of branches){const q=closest(x,y,z,b.id),first=b.segments[0],inlet=[x,y,z].reduce((s,v,k)=>s+(v-first.a[k])*first.tangent[k],0);tube=Math.max(tube,Math.min(q.clearance,inlet));}
   return Math.min(Math.max(main,tube),legacy.end-x)-dp/2;
  }
  function release(off,branch=0){
   const b=branches[branch],distance=cad&&p.releaseMode!=='distance'?b.length-.0005:Math.min(p.release/1000,b.length-.0005),s=b.length-distance,seg=b.segments.find(v=>s<=v.s+v.L)||b.segments.at(-1),t=s-seg.s;
   // A stable transverse frame for each curved passage. Offset never becomes
   // an axial displacement on a non-planar branch.
   let e=[-seg.tangent[1],seg.tangent[0],0],norm=Math.hypot(...e);if(norm<1e-5){e=[0,-seg.tangent[2],seg.tangent[1]];norm=Math.hypot(...e);}
   return seg.a.map((v,k)=>v+t*seg.tangent[k]+off*e[k]/norm);
  }
  return {...legacy,branches,closest,clearance,release,cad,source,dimension:3};
 }
 function makeField(g,p,f,legacyGeometry){
  const candidate=p.flowPhysics==='phase2',notional=JET.notional(f,g.dn/1000),jetSpeed=candidate&&notional.valid?notional.u:f.uJet,baseGeo=geometry(g,p,legacyGeometry),geo=p.wallModel==='cad'?WALL.attach(baseGeo,g.key):baseGeo,dn=geo.dn,b0=(candidate&&notional.valid?notional.diameter:dn)/Math.sqrt(2),bleedX=g.xh!=null?(g.xh-(g.xn||0))/1000:geo.end*.5;
  const smoothJunction=candidate&&geo.branches.every(b=>Math.abs(b.segments.at(-1).tangent[0])<.1);
  const fractions=geo.branches.map((b,i)=>f.circuit?.flowUsed?f.circuit.branches[i]?.fraction||0:1/geo.branches.length);
  const entries=geo.branches.map(b=>b.entry[0]-b.radius);
  // Linked diffusivity closure: changing turbulence intensity or mixing length
  // also changes spreading. The Mach reduction concerns the shear layer only.
  const Mc=Math.abs(jetSpeed-f.ms/(f.inletP/(R*f.Tin)*area(geo.db)))/(2*Math.sqrt(1.4*R*f.Tin));
  const spread=clamp(p.spread*(g.spreadScale??1)*(p.intensity/.18)*(p.lengthScale/.12)*(p.spreadingLaw==='compressible-fit'?.23+.77*Math.exp(-3.4*Mc*Mc):1/Math.sqrt(1+.25*Mc*Mc)),.02,.35),h=dn*.035;
  const referenceRadii=Array.from({length:101},(_,i)=>geo.radius(geo.end*i/100)),confinement=Object.freeze({...JET.confinement(f,g.d/2000),constantArea:Math.max(...referenceRadii)-Math.min(...referenceRadii)<1e-6});
  let noiseScale=1;
  const ingress=x=>fractions.reduce((s,w,i)=>s+w*smooth((x-entries[i]+geo.branches[i].radius)/(2*geo.branches[i].radius)),0);
  const axialCache=new Map();
  function axial(x){
   x=clamp(x,0,geo.end);if(axialCache.has(x))return axialCache.get(x);const r=geo.radius(x),b=b0+spread*x,F=ingress(x),H=smooth((x-bleedX+dn)/(2*dn)),mass=f.mp+f.ms*F+f.mh*H,powder=f.md*F;
   const totalT=(f.mp*f.T0+f.ms*F*f.Tin+f.mh*H*293.15)/mass,split=clamp(bleedX,dn,geo.end*.8);
   const middle=f.network?.bleedPressure||f.throatP,pressure=x<split?f.portP+(middle-f.portP)*x/split:middle+(f.po-middle)*(x-split)/(geo.end-split);
   const ec=(1+powder/mass)/(2*cp)*(mass*R/(pressure*Math.PI*r*r))**2,temp=2*totalT/(1+Math.sqrt(1+4*ec*totalT));
   const rho=pressure/(R*temp),mean=mass/(rho*Math.PI*r*r),fb=b*b/(r*r)*(-Math.expm1(-r*r/(b*b))),amp0=Math.min(jetSpeed*b0/b,mean/Math.max(fb,1e-9)*.98);
   const recirc=candidate?JET.recirculation({x,r,b,rho,mean,fb,amp:amp0,F,H},f,notional.valid?notional.diameter:dn,spread,confinement):{amp:amp0,reverse:0,Ct:null,demand:0},amp=recirc.amp;
   // Analytic Gaussian moments; normalization uses gas flux, not area alone.
   const moment=k=>b*b/(k*r*r)*(-Math.expm1(-k*r*r/(b*b))),e2=moment(2),e3=moment(3),variance=e2-fb*fb;
   const beta=1+(amp/mean)**2*variance,alpha=1+3*(amp/mean)**2*variance+(amp/mean)**3*(e3-3*fb*e2+2*fb**3);
   const result=Object.freeze({recirculation:Object.freeze(recirc),x,r,b,mass,powder,temp,pressure,rho,mean,fb,amp,turbulenceAmp:amp0,beta,alpha,fluxGaussian:fb+amp/mean*variance,F,H});
   if(axialCache.size>=2048)axialCache.delete(axialCache.keys().next().value);axialCache.set(x,result);return result;
  }
  const junctionLimit=Math.max(0,...geo.branches.map(b=>b.entry[0]+b.radius)),junctionScale=[],junctionPowderScale=[];
  function sectionScale(values,x){if(!smoothJunction||!values.length||x>=junctionLimit)return 1;const u=clamp(x/junctionLimit,0,1)*(values.length-1),i=Math.min(values.length-2,Math.floor(u));return values[i]+(values[i+1]-values[i])*(u-i);}
  const mainScale=x=>sectionScale(junctionScale,x);
  const integral=(a,y)=>a.rho*(a.mean*y*y/2+a.amp*(a.b*a.b/2*(-Math.expm1(-y*y/(a.b*a.b)))-a.fb*y*y/2));
  function sample(x,y,z=0,branchIndex=null){
   const a=axial(x),rad=Math.hypot(y,z),eta=rad/a.b,gauss=Math.exp(-eta*eta),ux0=a.mean+a.amp*(gauss-a.fb),am=axial(x-h),ap=axial(x+h),dx=Math.max(ap.x-am.x,1e-12);
   const ur=rad>1e-8?-(integral(ap,rad)-integral(am,rad))/(dx*a.rho*rad):0;
   const insideMain=x>=0&&x<=geo.end&&rad<a.r;
   const q=smoothJunction?geo.closest(x,y,z,branchIndex):insideMain?{branch:geo.branches[0],segment:geo.branches[0].segments[0],s:0,clearance:0}:geo.closest(x,y,z,branchIndex);
   // Keep suction-directed transport through the main-bore overlap. The old
   // hard insideMain switch turned particles sideways at the top of the bore,
   // before they could reach the jet. Join smoothly at the branch centreline end.
   const w=smoothJunction?smooth(q.clearance/(q.segment.radius*.12))*smooth((q.branch.length-q.s)/q.segment.radius):insideMain?0:clamp(q.clearance/(q.segment.radius*.12),0,1);
   const branchFlow=f.ms*fractions[q.branch.id],common=q.segment===q.branch.segments[0]&&geo.source.common.lengthMM>0,branchMass=common?f.ms:branchFlow;
   const inletPressure=f.inletP+(f.portP-f.inletP)*q.s/q.branch.length,branchRho=inletPressure/(R*f.Tin),us=branchMass/(branchRho*area(q.segment.radius*2));
   const contraction=candidate?us*(q.radiusSlope||0)/q.segment.radius:0,branchU=q.segment.tangent.map((v,i)=>us*v+contraction*(q.offset?.[i]||0));
   const ux=candidate?(1-w)*ux0*mainScale(x)+w*branchU[0]:(1-w)*ux0+w*us*q.segment.tangent[0],uy=candidate?(1-w)*(rad?ur*y/rad:0)+w*branchU[1]:(1-w)*(rad?ur*y/rad:0)+w*us*q.segment.tangent[1],uz=candidate?(1-w)*(rad?ur*z/rad:0)+w*branchU[2]:(1-w)*(rad?ur*z/rad:0)+w*us*q.segment.tangent[2];
   const development=Math.sqrt(1-Math.exp(-((Math.max(0,x)/((g.developmentD??.7)*dn))**2))),sigmaMain=(g.turbulenceScale??1)*p.intensity*a.turbulenceAmp*development*Math.exp(-(((eta-.7)/.65)**2));
   const sigma=noiseScale*(candidate?Math.hypot(Math.sqrt(1-w)*sigmaMain,Math.sqrt(w)*.15*p.intensity*us):(1-w)*sigmaMain+w*.15*p.intensity*us),ell=Math.max(dn*.015,(1-w)*p.lengthScale*Math.min(a.b,a.r)+w*.1*q.segment.radius*2),k=1.5*sigma*sigma,eps=.9*k**1.5/ell;
   const pressure=(1-w)*a.pressure+w*inletPressure,temp=(1-w)*a.temp+w*f.Tin,nu=mu(temp)*R*temp/pressure;
   const meanEntry=entries.reduce((s,v,i)=>s+v*fractions[i],0),contrast=.9*Math.exp(-Math.max(0,x-meanEntry)/(Math.max(1e-5,geo.db*2)));
   const theta=Math.atan2(z,y),azimuth=fractions.reduce((s,v,i)=>s+v*Math.cos(theta-geo.branches[i].angle),0);
   const scalarShape=a.recirculation.reverse>0?1:(1+contrast*(a.fluxGaussian-gauss))*(1+.7*contrast*azimuth);
   const mainConcentration=a.rho*a.powder/a.mass*scalarShape,branchConcentration=branchRho*f.md/f.ms,baseConcentration=(1-w)*mainConcentration+w*branchConcentration,lambda=sectionScale(junctionPowderScale,x)-1,concentration=candidate?baseConcentration*Math.exp(lambda*(1-w)):(1-w)*a.rho*a.powder/a.mass*scalarShape+w*branchRho*f.md/f.ms;
   return {ux,uy,uz,baseConcentration,mainConcentration,ux0,mainPowderDelta:mainConcentration*((1-w)*ux-mainScale(x)*ux0),branchPowderFlux:w*branchConcentration*ux,inletWeight:w,mainMassFlux:pressure/(R*temp)*(1-w)*ux0,branchMassFlux:pressure/(R*temp)*w*branchU[0],k,eps,G:Math.sqrt(eps/nu),meanGradient:Math.abs(2*a.amp*rad/a.b**2*gauss),sigma,ell,nu,p:pressure,T:temp,rho:pressure/(R*temp),powderConcentration:concentration,bulkPowderConcentration:a.rho*a.powder/a.mass,scalarShape,b:a.b,r:a.r,TL:NUM.lagrangianTime(sigma,eps)};
  }
  // Normalize the blended axial mass flux section by section. The branch
  // direction is retained; only the main-flow component is scaled. Uniform
  // area quadrature and linear interpolation are checked independently below.
  if(smoothJunction&&junctionLimit>0){
   const ns=Math.max(24,Math.ceil(junctionLimit/(dn*.04))),scales=[];
   for(let ix=0;ix<=ns;ix++){const x=junctionLimit*ix/ns,a=axial(x),nr=96,nt=128;let main=0,branch=0;
    for(let ir=0;ir<nr;ir++)for(let it=0;it<nt;it++){const r=a.r*Math.sqrt((ir+.5)/nr),theta=2*Math.PI*(it+.5)/nt,v=sample(x,r*Math.cos(theta),r*Math.sin(theta));main+=v.mainMassFlux-a.rho*(a.mean+a.amp*(Math.exp(-((r/a.b)**2))-a.fb));branch+=v.branchMassFlux;}
    const cell=Math.PI*a.r*a.r/(nr*nt),scale=(a.mass-branch*cell)/(a.mass+main*cell);if(!(scale>0&&Number.isFinite(scale)))throw Error('Junction flux normalization has no forward-flow solution.');scales.push(scale);
   }
   junctionScale.push(...scales);
   if(f.md>0){const cs=[];for(let ix=0;ix<=ns;ix++){
    const x=junctionLimit*ix/ns,a=axial(x),nr=96,nt=128,cell=Math.PI*a.r*a.r/(nr*nt),terms=[];if(!(a.powder>0)){cs.push(1);continue;}
    let baseline=0;for(let ir=0;ir<nr;ir++)for(let it=0;it<nt;it++){const r=a.r*Math.sqrt((ir+.5)/nr),theta=2*Math.PI*(it+.5)/nt,v=sample(x,r*Math.cos(theta),r*Math.sin(theta));terms.push([v.baseConcentration*v.ux*cell,1-v.inletWeight]);baseline+=v.mainConcentration*v.ux0*mainScale(x)*cell;}
    // Positive exponential scalar normalization preserves the undiluted branch
    // (w=1). This avoids a negative main-stream concentration at high loading.
    let lambda=0;const correction=a.powder*mainScale(x)-baseline;
    for(let j=0;j<16;j++){let flux=correction*Math.exp(lambda),derivative=flux;for(const [v,q] of terms){const f=v*Math.exp(lambda*q);flux+=f;derivative+=q*f;}if(Math.abs(flux/a.powder-1)<1e-10)break;lambda-=clamp((flux-a.powder)/derivative,-2,2);}
    if(!Number.isFinite(lambda))throw Error('Junction scalar normalization is non-finite.');cs.push(1+lambda);
   }junctionPowderScale.push(...cs);}


  }
  let lossPower=0;const nx=60,nr=24,dx=geo.end/nx;
  for(let i=0;i<nx;i++){const x=(i+.5)*dx,r=geo.radius(x),dr=r/nr;for(let j=0;j<nr;j++){const y=(j+.5)*dr,v=sample(x,-y);lossPower+=v.rho*v.eps*2*Math.PI*y*dr*dx;}}
  const out=axial(geo.end),kinetic=.5*(f.mg+f.md)*out.mean**2*out.alpha,available=Math.max(0,(candidate&&notional.valid?f.power-.5*f.mp*f.uJet**2+notional.power:f.power)-kinetic);
  if(lossPower>available&&lossPower>0)noiseScale=(available/lossPower)**(1/3);
  let peak=0,peakX=0,peakY=0,peakNu=mu(f.T)*R*f.T/f.po;
  for(let i=0;i<=180;i++){const x=dn*1.001+(geo.end-dn*1.001)*i/180,r=geo.radius(x);for(let j=0;j<=48;j++){const y=-r+2*r*j/48,v=sample(x,y);if(v.G>peak){peak=v.G;peakX=x;peakY=y;peakNu=v.nu;}}}
  return {smoothJunction,junctionPowderScale,junctionScale,junctionLimit,candidate,notional,confinement,geo,sample,axial,peak,peakX,peakY,peakNu,noiseScale,lossPower:lossPower*noiseScale**3,available,kinetic,spread,convectiveMach:Mc,branchFractions:fractions,modelVersion:2,scalarClosure:'Flux-normalized radial/azimuthal mixing reconstruction; passive powder, no inferred deposition.'};
 }
 return {geometry,makeField};
})();
if(typeof module!=='undefined')module.exports=EductorSpatial;

'use strict';
// Geometric encounter kernels, SI units. No sticking, growth or depletion here.
// Fuchs transition kernel: WRF/MOSAIC module_mosaic_coag.F, lines 894–930.
// WZW relative velocity: Ayala, Rosa & Wang (2008), equations 22–24;
// theta = tau / (0.4 u'^2 / epsilon), following Wang et al. (2000).
const ExposureCollisions=(()=>{
 const NUM=typeof ExposureNumerics!=='undefined'?ExposureNumerics:require('./numerics.js');
 const kB=1.380649e-23, Rgas=287.05;
 const SCREEN={radiusEta:.1,turbulentMach:.3,particleRe:.1,densityRatio:100,tracerSt:.1};
 const KEYS=['brownian','turbulent','inertial','unsupported','tracer'];
 function properties(d,rhop,p,T){
  if(![d,rhop,p,T].every(v=>Number.isFinite(v)&&v>0))throw Error('Collision kernels require positive finite diameter, density, absolute pressure and temperature.');
  const {mu,rho,lambda}=NUM.gas(p,T),Kn=2*lambda/d,Cc=1+Kn*(1.257+.4*Math.exp(-1.1/Kn)),mass=rhop*Math.PI*d**3/6;
  const D=kB*T*Cc/(3*Math.PI*mu*d),c=Math.sqrt(8*kB*T/(Math.PI*mass)),ell=8*D/(Math.PI*c),q=ell/d;
  // Algebraically rationalized Fuchs length. Avoid subtracting nearly equal
  // cubes for q << 1 or q >> 1; the small-q series avoids final cancellation.
  const h=Math.hypot(1,q),delta=d*(q<1e-4?q/2+q*q/3-q**3/8:2*((1+q)**2+(1+q)*h+h*h)/(3*(1+q+h))-1);
  return {d,mass,D,c,delta,Kn,Cc,lambda,mu,rho,nu:mu/rho,tau:rhop*d*d*Cc/(18*mu)};
 }
 function fuchs(a,b){const R=(a.d+b.d)/2,D=a.D+b.D,c=Math.hypot(a.c,b.c);return 4*Math.PI*R*D/(R/(R+Math.hypot(a.delta,b.delta))+4*D/(R*c));}
 function brownian(d1,d2,rhop,p,T){return fuchs(properties(d1,rhop,p,T),properties(d2,rhop,p,T));}
 function saffmanTurner(d1,d2,G){return Math.sqrt(8*Math.PI/15)*((d1+d2)/2)**3*Math.max(0,G);}
 function accelerationVariance(tau1,tau2,sigma,eps,nu,{TL:givenTL,convention='zhou2001'}={}){
  if(!(sigma>0&&eps>0&&nu>0&&tau1>0&&tau2>0))return 0;
  const TL=givenTL??NUM.lagrangianTime(sigma,eps),a=tau1/TL,b=tau2/TL,s=a+b,ab=a*b,phi=Math.max(a/b,b/a),gamma=(convention==='restatement2025'?1:phi)*.183*sigma*sigma/Math.sqrt(eps*nu),A=(1+a)*(1+b),B=(1+gamma*a)*(1+gamma*b);
  // WZW correction to the Kruis–Kusters acceleration variance. Rationalizing
  // gamma/(gamma-1)*(1/A-1/B) removes the removable singularity at gamma=1.
  const decorrelation=((a-b)**2+4*ab*(ab/(A*(1+Math.sqrt((1+s)/A)))))/s;
  return (convention==='restatement2025'?1/3:1)*(1+.6*Math.exp(-((phi-1)**1.5)))*sigma*sigma*decorrelation*gamma*(s+(gamma+1)*ab)/(A*B);
 }
 // Zhou, Wexler & Wang (2001), p.102. DNS fit at contact R/eta=1;
 // at other contact radii this is a sensitivity upper bound, not validation.
 // Below St=.1 the empirical tanh tails are tapered to recover g→1.
 // This asymptotic guard is our extension, not a fitted Zhou coefficient.
 function rdfMono(alpha,Rlambda){const z0=.5*(1+Math.tanh((alpha-.5)/.25)),z1=.5*(1+Math.tanh((alpha-1.25)/.1)),z2=.5*(1+Math.tanh((alpha-6.5)/2.5)),y0=18*alpha*alpha,y1=.36*alpha**2.5*Math.exp(-(alpha**2.5)),y2=.24*Math.exp(-.5*alpha),y3=.013*Math.exp(-.07*alpha);const excess=y0*(1-z0*z0)+Rlambda*z0*z0*(y1*(1-z1)+y2*z1+y3*z2);return 1+excess*Math.min(1,(alpha/.1)**2);}
 function rdf(tau1,tau2,G,sigma,nu){if(!(G>0&&sigma>0))return 1;const psi=Math.max(tau1/tau2,tau2/tau1),correlation=Math.min(1,2.6*Math.exp(-psi)+.205*Math.exp(-.0206*psi)*.5*(1+Math.tanh(psi-3))),Rlambda=Math.sqrt(15)*sigma*sigma/(nu*G);return 1+correlation*Math.sqrt(Math.max(0,rdfMono(tau1*G,Rlambda)-1)*Math.max(0,rdfMono(tau2*G,Rlambda)-1));}
 function abrahamson(R,tau1,tau2,sigma,TL){return Math.sqrt(8*Math.PI)*R*R*sigma*Math.sqrt(1/(1+tau1/TL)+1/(1+tau2/TL));}
 function pair(a,d1,d2,rhop,{legacy=false,tracerOnly=false,slip=0,convention='zhou2001',extended=false,clustering=false}={}){
  const i=properties(d1,rhop,a.p,a.T),j=d1===d2?i:properties(d2,rhop,a.p,a.T),G=Math.max(0,a.G||0),nu=i.nu,eps=G*G*nu,sigma=Math.max(0,a.sigma||0),R=(d1+d2)/2,eta=G>0?Math.sqrt(nu/G):Infinity;
  const st=Math.max(i.tau,j.tau)*G,size=R/eta,Mt=Math.sqrt(3)*sigma/Math.sqrt(1.4*Rgas*a.T),TL=a.TL??NUM.lagrangianTime(sigma,eps);
  const relativeSlip=q=>sigma>0?Math.sqrt(3)*sigma*Math.sqrt(q.tau/(TL+q.tau)):0;
  // Estimated turbulent slip for both species, plus measured tracked-particle
  // slip when available. The partner slip remains a local statistical estimate.
  const maxRe=Math.max(i.d*Math.max(Math.abs(slip),relativeSlip(i)),j.d*relativeSlip(j))/nu;
  const tracer=saffmanTurner(d1,d2,G),variance=legacy||tracerOnly?0:accelerationVariance(i.tau,j.tau,sigma,eps,nu,{TL,convention}),inertialKernel=Math.sqrt(8*Math.PI)*R*R*Math.sqrt(variance);
  // Combine turbulent shear and acceleration variances, not their kernels.
  const kinetic=!legacy&&!tracerOnly&&extended&&(st>100||size>1),g12=!legacy&&!tracerOnly&&clustering&&!kinetic?rdf(i.tau,j.tau,G,sigma,nu):1,unclustered=kinetic?abrahamson(R,i.tau,j.tau,sigma,TL):Math.hypot(tracer,inertialKernel),turbulent=unclustered*g12,brownianKernel=legacy?0:fuchs(i,j),rate=brownianKernel+turbulent;
  const flags=[];if(!kinetic&&size>SCREEN.radiusEta)flags.push('finite-size');if(Mt>SCREEN.turbulentMach)flags.push('compressible-turbulence');if(!kinetic&&maxRe>SCREEN.particleRe)flags.push('non-Stokes-drag');if(rhop/i.rho<SCREEN.densityRatio)flags.push('density-ratio');if(G>0&&!(sigma>0)&&!legacy&&!tracerOnly)flags.push('missing-turbulence-variance');if((legacy||tracerOnly)&&st>SCREEN.tracerSt)flags.push('tracer-inertia');
  if(kinetic&&Math.min(i.tau,j.tau)/TL<10)flags.push('kinetic-decorrelation-unverified');if(!legacy&&!tracerOnly&&clustering&&!kinetic&&Math.abs(size-1)>.1)flags.push('RDF-contact-radius-extrapolation');const unsupported=flags.length?turbulent:0;
  return {TL,convention,g12,kinetic,unclustered,rate,brownian:brownianKernel,turbulent,inertial:Math.max(0,turbulent-(kinetic?0:tracer)),tracer:kinetic?0:tracer,screened:brownianKernel+(flags.length?0:turbulent),unsupported,maxSt:st,maxSize:size,maxRe,Mt,flags,regime:st<=SCREEN.tracerSt?'tracer':'inertial'};
 }
 return {rdfMono,rdf,abrahamson,properties,fuchs,brownian,saffmanTurner,accelerationVariance,pair,SCREEN,KEYS,label:model=>model==='legacy'?'Legacy Saffman–Turner (no Brownian)':model==='tracer-brownian'?'Fuchs Brownian + Saffman–Turner':'Fuchs Brownian + WZW turbulence',version:'fuchs-zhou2001-st-reference-v2',conventions:{zhou2001:'Zhou et al. 2001 pp.101–102: gamma includes phi; radial acceleration variance has no extra 1/3.',restatement2025:'Sensitivity only: 2025 restatement without phi and with 1/3. Not the original Zhou convention.'}};
})();
if(typeof module!=='undefined')module.exports=ExposureCollisions;

'use strict';
/* Reduced particle-exposure model. All SI internally. No CFD fields or fitted slide results. */
const Exposure=(()=>{
 const dispersion=typeof ExposureDispersion!=='undefined'?ExposureDispersion:require('./dispersion.js');
 const numerics=typeof ExposureNumerics!=='undefined'?ExposureNumerics:require('./numerics.js');
 const network=typeof EductorNetwork!=='undefined'?EductorNetwork:typeof require==='function'?require('./network.js'):null;
 const upgrades=typeof EductorUpgrades!=='undefined'?EductorUpgrades:typeof require==='function'?require('../upgrades/physics.js'):null;
 const spatial=typeof EductorSpatial!=='undefined'?EductorSpatial:typeof require==='function'?require('../upgrades/spatial.js'):null;
 const collisions=typeof ExposureCollisions!=='undefined'?ExposureCollisions:typeof require==='function'?require('./collisions.js'):null;
 const population=typeof ExposurePopulation!=='undefined'?ExposurePopulation:require('./population.js');
 const doseKeys=collisions.KEYS,cumulativeKeys=['hazard','probability','screened','brownian','turbulent','unsupported'];
 const Rg=287.05,gamma=1.4,cp=gamma*Rg/(gamma-1),ATM=101325;
 const clamp=(x,a,b)=>Math.max(a,Math.min(b,x)),smooth=x=>(x=clamp(x,0,1),x*x*(3-2*x)),area=d=>Math.PI*d*d/4;
 const mu=T=>1.716e-5*(T/273.15)**1.5*(273.15+110.4)/(T+110.4);
 function particle(d,rhop,p,T,slip=0){
  const {mu:viscosity,rho,lambda}=numerics.gas(p,T),Kn=2*lambda/d;
  const Cc=1+Kn*(1.257+.4*Math.exp(-1.1/Kn)),Re=rho*d*Math.abs(slip)/viscosity,f=Re<1000?1+.15*Re**.687:.44*Re/24;
  return {tau:rhop*d*d*Cc/(18*viscosity*f),lambda,Cc,Re,Kn,rho,mu:viscosity,nu:viscosity/rho};
 }
 function nozzle(p0,p,T,d,Cd){
  if(!(p0>p&&p>0&&T>0))return null;
  const crit=(2/(gamma+1))**(gamma/(gamma-1)),choked=p/p0<=crit,pe=choked?p0*crit:p;
  const Te=T*(pe/p0)**((gamma-1)/gamma),ue=Math.sqrt(2*cp*(T-Te)),mp=Cd*area(d)*pe/(Rg*Te)*ue;
  return {mp,choked,uJet:Math.sqrt(2*cp*T*(1-(p/p0)**((gamma-1)/gamma))),power:mp*cp*T*(1-(p/p0)**((gamma-1)/gamma))};
 }
 function flow(g,p,legacy){
  p={...p,po:p.ps}; // One external reservoir pressure; internal stage pressures remain solved.
  if(![p.p0,p.ps,p.po,p.temp,p.feed].every(Number.isFinite)||p.ps<=0||p.po<=0||p.temp<=-273.15||p.feed<0)return {ok:false,reason:'Positive absolute pressures and temperature, and nonnegative feed, are required.'};
  const po=p.po*1e5,ps=p.ps*1e5,T=p.temp+273.15,p0=p.p0*1e5;
  if(p.mode==='prescribed'){
   const chamber=p.chamber*1e5,n=nozzle(p0,chamber,T,g.dn/1000,g.family==='british'?.94:.82);
   if(!n||!(p.air>0)||p.bleed<0)return {ok:false,reason:'Prescribed mode needs motive pressure above chamber pressure and positive feeder air.'};
   const ms=p.air/1000,mh=g.family==='british'?p.bleed/1000:0;
   return {ok:true,manual:true,caution:true,reason:'Prescribed air flows: the pressure / entrainment balance has not been solved.',...n,ms,mh,mg:n.mp+ms+mh,md:p.feed/3600,p0,ps,po,Tin:293.15,inletP:ps,portP:chamber,throatP:chamber,T0:T,T:(n.mp*T+(ms+mh)*293.15)/(n.mp+ms+mh)};
  }
  if(p.modelVersion===2){const f=upgrades.solve(g,p),bad=numerics.finiteScalars(f);return bad.length?{ok:false,constraints:['non-finite-output'],reason:'Non-finite flow output: '+bad.join(', ')}:f;}
  if(!(p0>po))return {ok:false,reason:'The existing solver requires motive stagnation pressure above outlet pressure.'};
  legacy.setBoundary(ps/ATM,po/ATM);let r,inletP=ps,trainNeeds=false;
  if(g.key.startsWith('nov'))r=legacy.solveCase(g,{suction:ps/ATM,outlet:po/ATM,pg:(p0-po)/1e5,feed:p.feed,temp:p.temp,loss:p.loss});
  else{
   const a=legacy.solveA(p0-po,p.feed,0,0,T),b=legacy.solveB(p0-po,p0-po,p.feed,0,0,T,0),c=legacy.solveC(p0-po,p.feed,0,T);
   const all=legacy.prepareResults(a,b,c);r=all[{revA:0,revB1:1,revB2:2,revC:3}[g.key]];if(g.key==='revB2')inletP=po+all[1].pt;
   trainNeeds=g.key.startsWith('revB')&&all.slice(1,3).some(v=>network.assess(v,{...g,db:v.st.db*1000,nports:1,df:v.st.db*1000},p).needed);
  }
  const gate=network.assess(r,g,p);
  if(gate.needed||trainNeeds){const f=network.solve(g,p);return {...f,gate:{needed:true,reason:trainNeeds?'At least one train stage needs a compressible estimate':gate.reason},raw:f.ok?network.legacyResult(f,g,r,po):r};}
  if(!r.supported)return {ok:false,reason:r.reason,raw:r};
  return {ok:true,manual:false,caution:r.caution,reason:r.caution?'Existing flow estimate is compressibility-sensitive.':'Existing reduced flow estimate passes its screening checks.',mp:r.mp,ms:r.ms,mh:r.mh||0,mg:r.mg,md:r.md,p0,ps,po:g.key==='revB1'?po+r.pe:po,systemOutletP:po,Tin:r.st?.Ts||293.15,inletP,portP:po+(r.pCh??r.p1),throatP:po+(r.pThroat??r.pMix??r.pe),T0:T,T:r.Tmix,uJet:r.n.uj,power:r.n.power,choked:r.n.choked,raw:r};
 }
 function geometry(g,release){
  const pts=g.profile.map(([x,r])=>[x/1000,r/1000]),end=g.end/1000,db=g.db/1000,dn=g.dn/1000;
  const angle=g.angle*Math.PI/180,sn=Math.sin(angle),cs=Math.cos(angle),R0=g.chamber/2000;
  function radius(x){x=clamp(x,0,end);let lo=0,hi=pts.length-1;while(hi-lo>1){const k=(lo+hi)>>1;if(pts[k][0]<x)lo=k;else hi=k;}const a=pts[lo],b=pts[hi],f=b[0]>a[0]?clamp((x-a[0])/(b[0]-a[0]),0,1):0;return a[1]+f*(b[1]-a[1]);}
  const cap=release/1000+.008;
  function clearance(x,y,dp=0){
   const h=1e-5,rr=radius(x),dr=(radius(x+h)-radius(x-h))/(2*h);
   const main=Math.min(x,(rr-Math.abs(y))/Math.sqrt(1+dr*dr),end-x)-dp/2;
   const s=-cs*x+sn*y,c=sn*x+cs*y;
   const branch=Math.min(s,cap-s,db/2-Math.abs(c))-dp/2;
   return Math.max(main,branch);
  }
  return {radius,end,db,dn,sn,cs,R0,cap,clearance,release:(off)=>[-cs*release/1000+sn*off,sn*release/1000+cs*off]};
 }
 function makeFieldCore(g,p,f){
  if(p.modelVersion===2)return spatial.makeField(upgrades.calibratedGeometry(g,p),p,f,geometry(g,p.release));
  const geo=geometry(g,p.release),b0=geo.dn/Math.sqrt(2),h=geo.dn*.035,bleedX=g.xh!=null?(g.xh-(g.xn||0))/1000:geo.end*.5;
  let noiseScale=1;
  function axial(x){
   x=clamp(x,0,geo.end);const r=geo.radius(x),b=b0+p.spread*x;
   const mixH=smooth((x-bleedX+geo.dn)/(2*geo.dn)),mass=f.mp+f.ms+f.mh*mixH;
   const totalT=(f.mp*f.T0+f.ms*f.Tin+f.mh*mixH*293.15)/mass;
   // Port -> mixing bore -> outlet interpolation; a reduced pressure field, not a nozzle/shock solution.
   const split=Math.min(g.L/1000,geo.end*.8),pressure=x<split?f.portP+(f.throatP-f.portP)*x/split:f.throatP+(f.po-f.throatP)*(x-split)/(geo.end-split);
   // Static temperature from adiabatic bulk energy including powder inertia.
   const energyC=(1+f.md/mass)/(2*cp)*(mass*Rg/(pressure*Math.PI*r*r))**2;
   const temp=f.compressible?2*totalT/(1+Math.sqrt(1+4*energyC*totalT)):totalT;
   const rho=pressure/(Rg*temp),mean=mass/(rho*Math.PI*r*r),fb=b*b/(r*r)*(-Math.expm1(-r*r/(b*b)));
   const amp=Math.min(f.uJet*b0/b,mean/Math.max(fb,1e-9));
   return {x,r,b,mass,temp,pressure,rho,mean,fb,amp};
  }
  const integral=(a,y)=>a.rho*(a.mean*y*y/2+a.amp*(a.b*a.b/2*(-Math.expm1(-y*y/(a.b*a.b)))-a.fb*y*y/2));
  function sample(x,y){
   const a=axial(x),eta=Math.abs(y)/a.b,uf=a.mean+a.amp*(Math.exp(-eta*eta)-a.fb);
   const am=axial(x-h),ap=axial(x+h),dx=Math.max(ap.x-am.x,1e-12);
   const radial=Math.abs(y)>1e-8?-(integral(ap,Math.abs(y))-integral(am,Math.abs(y)))/(dx*a.rho*Math.abs(y))*Math.sign(y):0;
   const s=-geo.cs*x+geo.sn*y,c=geo.sn*x+geo.cs*y,w=Math.abs(c)<geo.db/2&&s>0?smooth((s/geo.R0-.08)/.95):0;
   const feederRho=f.inletP/(Rg*f.Tin),us=f.ms/(g.nports*feederRho*area(geo.db));
   const ux=(1-w)*uf+w*geo.cs*us,uy=(1-w)*radial-w*geo.sn*us;
   const dev=1-Math.exp(-((Math.max(0,x)/(2*geo.dn))**2)),sigmaMain=p.intensity*a.amp*Math.sqrt(dev)*Math.exp(-(((eta-.7)/.65)**2));
   const sigma=noiseScale*((1-w)*sigmaMain+w*.15*p.intensity*us),ell=Math.max(geo.dn*.015,(1-w)*p.lengthScale*Math.min(a.b,a.r)+w*.1*geo.db);
   const k=1.5*sigma*sigma,eps=.9*k**1.5/ell;
   const pressure=(1-w)*a.pressure+w*f.inletP,temp=(1-w)*a.temp+w*f.Tin,nu=mu(temp)*Rg*temp/pressure;
   return {ux,uy,k,eps,G:Math.sqrt(eps/nu),sigma,ell,nu,p:pressure,T:temp,rho:pressure/(Rg*temp),powderConcentration:pressure/(Rg*temp)*f.md*((1-w)/a.mass+w/f.ms),b:a.b,r:a.r,TL:numerics.lagrangianTime(sigma,eps)};
  }
  // Bound the modeled turbulent dissipation by available motive mechanical power.
  let lossPower=0;const nx=60,nr=18,dx=geo.end/nx;
  for(let i=0;i<nx;i++){const x=(i+.5)*dx,r=geo.radius(x),dr=r/nr;for(let j=0;j<nr;j++){const y=(j+.5)*dr,v=sample(x,-y);lossPower+=v.rho*v.eps*2*Math.PI*y*dr*dx;}}
  const out=axial(geo.end),kinetic=.5*(f.mg+f.md)*out.mean*out.mean,available=Math.max(0,f.power-kinetic);
  if(lossPower>available&&lossPower>0)noiseScale=(available/lossPower)**(1/3);
  let peak=0,peakX=0,peakY=0,peakNu=mu(f.T)*Rg*f.T/f.po;
  for(let i=0;i<=180;i++){const x=geo.dn*1.001+(geo.end-geo.dn*1.001)*i/180;if(x<0)continue;const r=geo.radius(x);for(let j=0;j<=48;j++){const y=-r+(2*r)*j/48,v=sample(x,y);if(v.G>peak){peak=v.G;peakX=x;peakY=y;peakNu=v.nu;}}}
  return {geo,sample,axial,peak,peakX,peakY,peakNu,noiseScale,lossPower:lossPower*noiseScale**3,available,kinetic};
 }
 function makeField(g,p,f){const field=makeFieldCore(g,p,f),bad=numerics.finiteScalars(field);if(bad.length)throw Error('Non-finite field output: '+bad.join(', '));return field;}
 function random(seed){let x=seed>>>0;const u=()=>{x+=0x6D2B79F5;let t=x;t=Math.imul(t^t>>>15,t|1);t^=t+Math.imul(t^t>>>7,t|61);return ((t^t>>>14)>>>0)/4294967296;};return ()=>Math.sqrt(-2*Math.log(Math.max(u(),1e-12)))*Math.cos(2*Math.PI*u());}
 // Shared physical sections for dilution and particle shear. The first is the
 // release section in the suction branch; downstream sections are axial planes.
 function stationSections(field){
  const peak=field.peakX,end=field.geo.end;
  return [['Suction inlet',null],['Jet entry',0],['Peak shear',peak],['Mixing 7.5%',peak+.075*(end-peak),true],['Mixing 15%',peak+.15*(end-peak),true],['Mixing 22.5%',peak+.225*(end-peak),true],['Mixing 30%',peak+.3*(end-peak),true],['Midway',peak+.5*(end-peak)],['Mixing 70%',peak+.7*(end-peak),true],['Pre-outlet',peak+.85*(end-peak)],['Mixing 95%',peak+.95*(end-peak),true],['Outlet',end]].map(([label,x,diagnosticOnly],i)=>({id:i+1,label,x,...(diagnosticOnly?{diagnosticOnly:true}:{})}));
 }
 function sampleStationCrossings(field,d,sections,samples,a,b,dose=null){
  if(!(b.x>a.x))return;
  for(let i=1;i<sections.length-1;i++){
   const x=sections[i].x;if(samples[i]||x<a.x||x>b.x)continue;
   const w=(x-a.x)/(b.x-a.x),y=a.y+w*(b.y-a.y),z=(a.z||0)+w*((b.z||0)-(a.z||0));
   // Do not count a crossing of the same x-plane high in the suction branch.
   if(Math.hypot(y,z)>field.geo.radius(x)-d/2)continue;
   samples[i]={t:a.t+w*(b.t-a.t),x,y,z,G:field.sample(x,y,z).G,...(sections[i].diagnosticOnly?{diagnosticOnly:true}:{})};
   if(dose){
    const hazard=dose.hazard+integrateRate(dose.rate0,dose.rate1,b.t-a.t,w),screened=dose.screened+integrateRate(dose.screened0,dose.screened1,b.t-a.t,w);
    Object.assign(samples[i],{hazard,screened,probability:coagProbability(hazard,dose.alpha)});
    if(dose.components)for(const key of doseKeys)samples[i][key]=dose.components[key]+integrateRate(dose.start[key],dose.end[key],b.t-a.t,w);
   }
  }
 }
 function crossing(g0,g1,dt,threshold){if(Math.max(g0,g1)<threshold)return 0;if(Math.min(g0,g1)>=threshold)return dt;const f=(threshold-g0)/(g1-g0);return dt*(g1>g0?1-f:f);}
 function turnover(t,q,Gref,nu,d){const eta=Math.sqrt(nu/(q*Gref));return q*Gref*t*Math.min(1,(eta/d)**(2/3));}
 function track(field,p,d,offset,Gref,seed,noisy,keep=false,branchIndex=0){
  const geo=field.geo,rand=random(seed);let [x,y,z=0]=geo.release(offset,branchIndex),a=field.sample(x,y,z,branchIndex),vx=a.ux,vy=a.uy,vz=a.uz||0,t=0,wx=0,wy=0,wz=0;
  if(geo.clearance(x,y,d,z)<=0)throw Error('Release intersects the reduced passage wall; reduce the release offset.');
  const bins=241,horizon=p.horizon/1000,series=new Float64Array(bins),points=keep?[geo.dimension===3?[0,x,y,a.G,z]:[0,x,y,a.G]]:null;
  const sections=stationSections(field),stationSamples=sections.map(()=>null);stationSamples[0]={t:0,x,y,z,G:a.G,hazard:0,screened:0,probability:0,...Object.fromEntries(doseKeys.map(k=>[k,0]))};
  const hazardSeries=new Float64Array(bins),screenedSeries=new Float64Array(bins);let coag=collisionState(a,p,d),hazard=0,screened=0,shearDose=0,epsilonDose=0,maxSt=coag.maxSt,maxSize=coag.maxSize,maxMt=coag.Mt,maxConcentration=coag.concentration,maxVolume=coag.volumeFraction;
  const components=Object.fromEntries(doseKeys.map(k=>[k,0])),componentSeries=Object.fromEntries(doseKeys.map(k=>[k,new Float64Array(bins)]));let maxCollisionRe=coag.maxCollisionRe;
  const populationPath=p.populationBalance?[{t:0,a}]:null;let maxMach=0;
  let nextBin=1,half=0,quarter=0,maxRe=0,peak=0,steps=0,fate='censored',numericalReason=null,tSample=0,maxDt=0,minDt=Infinity;
  let normalized=[0,0,0];
  if(noisy&&p.dispersion==='well-mixed'&&p.initialFluctuation!=='zero'){normalized=normalized.map((_,i)=>numerics.normal(seed,i,0)*(p.antitheticSign??1));[wx,wy,wz]=normalized.map(q=>q*a.sigma);}
  const fluctuationScale=Math.max(a.sigma,Number.isFinite(field.peakX)?field.sample(field.peakX,field.peakY).sigma:0,1e-12);
  const bridge=p.commonNoise?numerics.brownian(seed,horizon,p.antitheticSign??1):null;
  series[0]=a.G;
  while(t<horizon&&steps<(p.maxSteps??2000000)){
   const speed=Math.max(Math.hypot(vx,vy,vz),Math.hypot(a.ux,a.uy,a.uz||0),1);
   const clearance=geo.clearance(x,y,d,z);
   if(clearance<=1e-9){const q=geo.dimension===3?geo.closest(x,y,z):null;fate=x>=geo.end-d/2-2e-6?'outlet':q&&q.s<=d/2+2e-6?'inlet':'wall';break;}
   let dt=numerics.timestep({maximum:p.dt*1e-6,Gref,G:a.G,b:a.b,r:a.r,db:geo.db,speed,remaining:horizon-t,TL:a.TL,noisy,clearance,adaptive:p.adaptiveDt===true});
   const grad=noisy&&p.dispersion==='well-mixed'?dispersion.gradients(field,x,y,z,a):null;
   if(grad&&p.inhomogeneousDt===true)dt=Math.min(dt,dispersion.stepLimit([wx,wy,wz],a,[vx-a.ux,vy-a.uy,vz-(a.uz||0)],grad));
   if(!(dt>1e-15)){fate='numerical';break;}
   maxDt=Math.max(maxDt,dt);minDt=Math.min(minDt,dt);
   const gaussian=axis=>bridge?bridge.increment(t,t+dt,axis)/Math.sqrt(dt):(p.antitheticSign??1)*rand();
   const ex=noisy?Math.exp(-dt/Math.max(a.TL,1e-12)):0,ns=noisy?a.sigma*Math.sqrt(-Math.expm1(-2*dt/Math.max(a.TL,1e-12))):0;
   const normalizedNext=grad?dispersion.advanceNormalized(normalized,a,[vx-a.ux,vy-a.uy,vz-(a.uz||0)],dt,gaussian,grad,{isotropic:p.referenceProtocol===true}):null;
   const mixed=normalizedNext?normalizedNext.map(q=>q*a.sigma):null;
   const nxw=mixed?mixed[0]:noisy?wx*ex+ns*gaussian(0):0,nyw=mixed?mixed[1]:noisy?wy*ex+ns*gaussian(1):0,nzw=mixed?mixed[2]:noisy&&geo.dimension===3?wz*ex+ns*gaussian(2):0,uz=(a.uz||0)+(wz+nzw)/2,ux=a.ux+(wx+nxw)/2,uy=a.uy+(wy+nyw)/2;
   // Reject catastrophic stochastic updates; never count their artificial
   // endpoint as a physical wall contact or exit. No velocity clipping.
   if(mixed&&(!normalizedNext.every(Number.isFinite)||Math.hypot(...mixed)>50*fluctuationScale)){fate='numerical';numericalReason='unresolved-dispersion-gradient';break;}
   const slip=Math.hypot(ux-vx,uy-vy,uz-vz),props=particle(d,p.density,a.p,a.T,slip),tau=props.tau;maxMach=Math.max(maxMach,slip/Math.sqrt(1.4*Rg*a.T));
   maxRe=Math.max(maxRe,props.Re);const ey=uy-9.81*tau,e=Math.exp(-dt/tau),om=-Math.expm1(-dt/tau);
   let xn=x+ux*dt+(vx-ux)*tau*om,yn=y+ey*dt+(vy-ey)*tau*om,zn=z+uz*dt+(vz-uz)*tau*om;
   let frac=1,hit=false;
   const contact=numerics.firstContact((xx,yy,zz)=>geo.clearance(xx,yy,d,zz),[x,y,z],[xn,yn,zn]);if(contact!==null){frac=contact;hit=true;}
   xn=x+frac*(xn-x);yn=y+frac*(yn-y);zn=z+frac*(zn-z);const tn=t+frac*dt,b=field.sample(xn,yn,zn,branchIndex),cb=collisionState(b,p,d,Math.hypot(b.ux-(ux+(vx-ux)*e),b.uy-(ey+(vy-ey)*e),(b.uz||0)-(uz+(vz-uz)*e))),elapsed=tn-t;
   sampleStationCrossings(field,d,sections,stationSamples,{t,x,y,z},{t:tn,x:xn,y:yn,z:zn},{hazard,screened,rate0:coag.rate,rate1:cb.rate,screened0:coag.screened,screened1:cb.screened,alpha:p.alpha??1,components,start:coag,end:cb});
   const nextHazard=hazard+integrateRate(coag.rate,cb.rate,elapsed),nextScreened=screened+integrateRate(coag.screened,cb.screened,elapsed);shearDose+=.5*(a.G+b.G)*elapsed;epsilonDose+=.5*(a.eps+b.eps)*elapsed;maxSt=Math.max(maxSt,cb.maxSt);maxSize=Math.max(maxSize,cb.maxSize);maxMt=Math.max(maxMt,cb.Mt);maxConcentration=Math.max(maxConcentration,cb.concentration);maxVolume=Math.max(maxVolume,cb.volumeFraction);maxCollisionRe=Math.max(maxCollisionRe,cb.maxCollisionRe);
   half+=crossing(a.G,b.G,tn-t,Gref/2);quarter+=crossing(a.G,b.G,tn-t,Gref/4);peak=Math.max(peak,a.G,b.G);
   while(nextBin<bins&&nextBin*horizon/(bins-1)<=tn){const tb=nextBin*horizon/(bins-1);const fracBin=(tb-t)/Math.max(tn-t,1e-15);series[nextBin]=a.G+(b.G-a.G)*fracBin;hazardSeries[nextBin]=hazard+integrateRate(coag.rate,cb.rate,elapsed,fracBin);screenedSeries[nextBin]=screened+integrateRate(coag.screened,cb.screened,elapsed,fracBin);for(const key of doseKeys)componentSeries[key][nextBin]=components[key]+integrateRate(coag[key],cb[key],elapsed,fracBin);nextBin++;}
   if(keep&&(tn>=tSample||hit)){points.push(geo.dimension===3?[tn,xn,yn,b.G,zn]:[tn,xn,yn,b.G]);tSample=tn+Math.max(horizon/1100,2e-6);}
   for(const key of doseKeys)components[key]+=integrateRate(coag[key],cb[key],elapsed);
   const populationPrevious=populationPath?{t,a}:null;
   vx=ux+(vx-ux)*Math.exp(-elapsed/tau);vy=ey+(vy-ey)*Math.exp(-elapsed/tau);vz=uz+(vz-uz)*Math.exp(-elapsed/tau);if(normalizedNext)normalized=normalizedNext;wx=normalizedNext?normalized[0]*b.sigma:nxw;wy=normalizedNext?normalized[1]*b.sigma:nyw;wz=normalizedNext?normalized[2]*b.sigma:nzw;x=xn;y=yn;z=zn;t=tn;a=b;coag=cb;hazard=nextHazard;screened=nextScreened;steps++;
   if(populationPath){
    const last=populationPath.at(-1),keys=['powderConcentration','G','eps','sigma','TL','p','T'],changed=(u,v)=>keys.some(k=>Math.abs((u[k]||0)-(v[k]||0))>(p.populationTolerance??.02)*Math.max(Math.abs(u[k]||0),Math.abs(v[k]||0),1e-20));
    if(hit||t-last.t>=(p.populationStep??20)*1e-6||changed(a,last.a)){
     // Preserve both sides of rapid field transitions; coarse uniform sampling
     // biases nonlinear kernels even when parcel integration conserves mass.
     if(changed(a,populationPrevious.a)&&populationPrevious.t>last.t)populationPath.push(populationPrevious);
     populationPath.push({t,a});
    }
   }
   if(hit){const q=geo.dimension===3?geo.closest(x,y,z):null,s=q?q.branch.length-q.s:-geo.cs*x+geo.sn*y;if(x>=geo.end-d/2-2e-6&&Math.hypot(y,z)<geo.radius(geo.end))fate='outlet';else if(q?q.s<=d/2+2e-6:s>=geo.cap-d/2-2e-6)fate='inlet';else fate='wall';break;}
  }
  if(steps>=(p.maxSteps??2000000)&&t<horizon)fate='numerical';
  if(keep&&points.at(-1)[0]!==t)points.push(geo.dimension===3?[t,x,y,a.G,z]:[t,x,y,a.G]);
  // A finite-size particle centre stops just inside the outlet plane.
  if(fate==='outlet')stationSamples[sections.length-1]={t,x,y,z,G:a.G,hazard,screened,probability:coagProbability(hazard,p.alpha??1),...components};
  hazardSeries.fill(hazard,nextBin);screenedSeries.fill(screened,nextBin);for(const key of doseKeys)componentSeries[key].fill(components[key],nextBin);
  if(populationPath&&populationPath.at(-1).t<t)populationPath.push({t,a});
  const growth=populationPath?evolvePopulation(populationPath,p,d,stationSamples):null;
  return {growth,maxMach,t,half,quarter,peak,maxRe,steps,maxDt,minDt:Number.isFinite(minDt)?minDt:0,meanDt:steps?t/steps:0,fate,numericalReason,contact:fate==='wall',deposited:null,branchIndex,dimension:geo.dimension||2,series,points,stationSamples,hazard,screened,probability:coagProbability(hazard,p.alpha??1),shearDose,epsilonDose,maxSt,maxSize,maxMt,maxConcentration,maxVolume,maxCollisionRe,...components,componentSeries,hazardSeries,screenedSeries};
 }
 function collisionOptions(p,slip=0){return {legacy:p.collisionModel==='legacy',tracerOnly:p.collisionModel==='tracer-brownian',slip,convention:p.kernelConvention||'zhou2001',extended:p.extendedKernels,clustering:p.clustering};}
 function evolvePopulation(path,p,d,stations=[]){
  const mixture=p.background==='mixture',ds=mixture?[p.small*1e-6,p.large*1e-6]:[d],shares=mixture?[p.smallShare??.5,1-(p.smallShare??.5)]:[1],pop=population.create(ds,shares,p.density,{bins:36});
  const stationGrowth=stations.map(()=>null);if(stations[0])stationGrowth[0]=population.summary(pop);const diameters=pop.m.map(m=>Math.cbrt(6*m/(p.density*Math.PI)));
  for(let j=1;j<path.length;j++){const start=path[j-1].t,end=path[j].t,breaks=[...stations.map((s,i)=>s&&!s.diagnosticOnly&&s.t>start&&s.t<=end?{t:s.t,index:i}:null).filter(Boolean),{t:end,index:null}].sort((a,b)=>a.t-b.t);let previous=start;for(const cut of breaks){
   for(let i=0;i<stations.length;i++){const station=stations[i];if(!station?.diagnosticOnly||station.t<=previous||station.t>cut.t)continue;const probe=structuredClone(pop),a=path[j-1].a,b=path[j].a,q={...a},fraction=(.5*(previous+station.t)-start)/Math.max(end-start,1e-30);for(const key of ['p','T','G','eps','sigma','TL','nu','rho','powderConcentration'])q[key]=a[key]+fraction*(b[key]-a[key]);population.advance(probe,Math.max(0,q.powderConcentration),station.t-previous,(i,k)=>collisions.pair(q,diameters[i],diameters[k],p.density,collisionOptions(p)).rate,p.alpha??1);stationGrowth[i]=population.summary(probe);}
   const a=path[j-1].a,b=path[j].a,q={...a},fraction=(.5*(previous+cut.t)-start)/Math.max(end-start,1e-30);for(const key of ['p','T','G','eps','sigma','TL','nu','rho','powderConcentration'])q[key]=a[key]+fraction*(b[key]-a[key]);population.advance(pop,Math.max(0,q.powderConcentration),cut.t-previous,(i,k)=>collisions.pair(q,diameters[i],diameters[k],p.density,collisionOptions(p)).rate,p.alpha??1);if(cut.index!==null)stationGrowth[cut.index]=population.summary(pop);previous=cut.t;}}
  return {...population.summary(pop),stations:stationGrowth,scope:'Sectional growth along a fixed transport parcel; growing sizes do not change the carrier field or parcel path.'};
 }
 function quantile(sorted,q){if(!sorted.length)return null;const x=(sorted.length-1)*q,i=Math.floor(x);return sorted[i]+(sorted[Math.min(i+1,sorted.length-1)]-sorted[i])*(x-i);}
 function stats(a){if(!a.length)return {n:0,mean:null,sd:null,se:null,ci:null,p05:null,p25:null,median:null,p75:null,p95:null,p99:null,p10:null,p90:null,cv:null,min:null,max:null};const n=a.length,mean=a.reduce((s,x)=>s+x,0)/n,sd=n>1?Math.sqrt(a.reduce((s,x)=>s+(x-mean)**2,0)/(n-1)):0,sorted=[...a].sort((a,b)=>a-b),se=sd/Math.sqrt(n),df=n-1,small=[0,12.7062,4.3027,3.1824,2.7764,2.5706,2.4469,2.3646,2.306,2.2622,2.2281],t=df<=10?small[df]:1.959964+2.372272/df+2.8225/(df*df)+2.556/(df**3);return {n,mean,sd,se,ci:n>1?[mean-t*se,mean+t*se]:null,p05:quantile(sorted,.05),p25:quantile(sorted,.25),median:quantile(sorted,.5),p75:quantile(sorted,.75),p95:quantile(sorted,.95),p99:quantile(sorted,.99),p10:quantile(sorted,.1),p90:quantile(sorted,.9),cv:mean!==0?sd/Math.abs(mean):null,min:sorted[0],max:sorted.at(-1)};}
 function weightedStats(values,weights){
  if(!values.length)return stats([]);
  const total=weights.reduce((a,b)=>a+b,0),w=weights.map(v=>v/total),mean=values.reduce((a,v,i)=>a+v*w[i],0),w2=w.reduce((a,b)=>a+b*b,0),effectiveN=1/w2,sd=Math.sqrt(values.reduce((a,v,i)=>a+w[i]*(v-mean)**2,0)/Math.max(1e-12,1-w2)),se=sd/Math.sqrt(effectiveN),sorted=values.map((v,i)=>({v,w:w[i]})).sort((a,b)=>a.v-b.v);
  const q=x=>{let c=0;for(const a of sorted){c+=a.w;if(c>=x)return a.v;}return sorted.at(-1).v;};
  const df=effectiveN-1,t=df>=1?1.959964+2.372272/df+2.8225/df**2+2.556/df**3:0;
  return {n:values.length,effectiveN,mean,sd,se,ci:effectiveN>1?[mean-t*se,mean+t*se]:null,p05:q(.05),p25:q(.25),median:q(.5),p75:q(.75),p95:q(.95),p99:q(.99),p10:q(.1),p90:q(.9),cv:mean?sd/Math.abs(mean):null,min:sorted[0].v,max:sorted.at(-1).v};
 }
 // A frozen-background encounter model. Full dose retains explicitly flagged
 // extrapolation; screened dose is a subtotal, never a corrected total.
 const collisionKernel=collisions.saffmanTurner;
 function collisionState(a,p,d,slip=0){
  const concentration=Math.max(0,a.powderConcentration||0),share=clamp(p.smallShare??.5,0,1),background=p.background==='mixture'?[[p.small*1e-6,share],[p.large*1e-6,1-share]]:[[d,1]],legacy=p.collisionModel==='legacy';
  const out={rate:0,screened:0,...Object.fromEntries(doseKeys.map(k=>[k,0])),maxSt:0,maxSize:0,maxCollisionRe:0,Mt:0,concentration,volumeFraction:concentration/p.density};
  for(const [dj,weight] of background){if(!weight)continue;const n=concentration*weight/(p.density*Math.PI*dj**3/6),q=collisions.pair(a,d,dj,p.density,collisionOptions(p,slip));
   for(const key of ['rate','screened',...doseKeys])out[key]+=q[key]*n;
   out.maxSt=Math.max(out.maxSt,q.maxSt);out.maxSize=Math.max(out.maxSize,q.maxSize);out.maxCollisionRe=Math.max(out.maxCollisionRe,q.maxRe);out.Mt=Math.max(out.Mt,q.Mt);
  }
  return out;
 }
 const coagProbability=(hazard,alpha=1)=>-Math.expm1(-clamp(alpha,0,1)*Math.max(0,hazard));
 const integrateRate=(a,b,dt,f=1)=>dt*(a*f+(b-a)*f*f/2);

 function wilson(k,n){const z=1.95996398454,d=1+z*z/n,c=(k/n+z*z/(2*n))/d,h=z*Math.sqrt((k/n)*(1-k/n)/n+z*z/(4*n*n))/d;return [Math.max(0,c-h),Math.min(1,c+h)];}
 function* simulation(g,p,f,options={}){
  const field=options.field||makeField(g,p,f),Gref=p.reference==='slide'?3829961.607121514:field.peak;
  if(!(Gref>0&&Number.isFinite(Gref)))return {ok:false,reason:'No positive turbulence reference: increase the turbulence intensity or provide an energetically feasible flow.'};
  const rows=[],offsets=[-1,0,1].map(v=>v*p.offset*g.db/2000),sizes=[p.small*1e-6,p.large*1e-6];
  for(let species=0;species<2;species++)for(let o=0;o<3;o++){
   if(options.groups&&!options.groups.includes(species*3+o))continue;
   const d=sizes[species],offset=offsets[o],adaptive=p.noise&&p.sampling==='precision',limit=p.noise?(adaptive?(p.maxHistories??4096):p.n):(field.geo.branches?.length||1),runs=[],mean=new Float64Array(241),sq=new Float64Array(241);let first;
   const fractions=field.branchFractions||[1],counts=fractions.map(()=>0),paired=p.noise&&p.antithetic===true;
   const seedFor=j=>(p.seed+species*103309+o*1009+(paired?Math.floor(j/2):j)*7919)>>>0;
   for(let j=0;j<limit;j++){
    const position=(((paired?Math.floor(j/2):j)+.5)*.618033988749895)%1;let c=0;const branch=!p.noise?j:(j<(paired?2:1)*fractions.length?Math.floor(j/(paired?2:1)):fractions.findIndex((v,i)=>(c+=v)>=position||i===fractions.length-1));
    const r=track(field,{...p,antitheticSign:paired&&j%2?-1:1},d,offset,Gref,seedFor(j),p.noise,counts[branch]===0,branch);counts[branch]++;runs.push(r);if(j===0)first=r;
    if((j+1)%(p.populationBalance||p.wallModel==='cad'?1:16)===0)yield {completed:(species*3+o)*limit+j+1,total:6*limit};
    if(adaptive&&(j+1)%128===0&&counts.every((n,i)=>n>0||fractions[i]===0)){
     runs.forEach(r=>r.weight=fractions[r.branchIndex]/counts[r.branchIndex]);
     const check=numerics.precision(runs,p.targetHalfWidth??.02);if(paired){check.halfWidth=(numerics.wilson(runs.reduce((s,r)=>s+(r.contact?r.weight:0),0),check.effectiveN/2)[1]-numerics.wilson(runs.reduce((s,r)=>s+(r.contact?r.weight:0),0),check.effectiveN/2)[0])/2;}
     if(check.halfWidth<=(p.targetHalfWidth??.02))break;
    }
   }
   const n=runs.length;for(const r of runs){r.weight=fractions[r.branchIndex]/counts[r.branchIndex];for(let b=0;b<241;b++){mean[b]+=r.series[b]*r.weight;sq[b]+=r.series[b]**2*r.weight;}}
   const summary=(items,fn)=>p.modelVersion===2?weightedStats(items.map(fn),items.map(r=>r.weight)):stats(items.map(fn));
   const sd=Array.from(mean,(v,i)=>n>1?Math.sqrt(Math.max(0,(sq[i]-v*v)/(p.modelVersion===2?1-runs.reduce((s,r)=>s+r.weight*r.weight,0):(n-1)/n))):0),wall=runs.filter(r=>r.fate==='wall').length;
   const sampleKeys=['t','half','quarter','peak','maxRe','maxDt','meanDt','maxMach','fate','hazard','screened','probability','shearDose','epsilonDose','maxSt','maxSize','maxMt','maxConcentration','maxVolume','maxCollisionRe',...doseKeys];
   const samples=runs.map((r,j)=>Object.fromEntries([['index',j+1],['seed',seedFor(j)],['antitheticSign',paired&&j%2?-1:1],['stationSamples',r.stationSamples],['branch',r.branchIndex],['powderFluxWeight',r.weight],['contact',r.contact],['deposited',r.deposited],['growth',r.growth],['numericalReason',r.numericalReason],...sampleKeys.map(k=>[k,r[k]])]));
   const stationShear=stationSections(field).map((section,i)=>({...section,total:n,...summary(runs.filter(r=>r.stationSamples[i]),r=>r.stationSamples[i].G)}));
   const stationCoagulation=Object.fromEntries(cumulativeKeys.map(key=>[key,stationSections(field).map((section,i)=>({...section,total:n,...summary(runs.filter(r=>r.stationSamples[i]),r=>r.stationSamples[i][key])}))]));
   const distribution=Object.fromEntries(sampleKeys.filter(k=>k!=='fate').map(k=>[k,summary(runs,r=>r[k])]));distribution.outletTime=summary(runs.filter(r=>r.fate==='outlet'),r=>r.t);
   const residence=numerics.residence(samples),precision=numerics.precision(runs,p.targetHalfWidth??.02);if(paired){const ci=numerics.wilson(runs.reduce((s,r)=>s+(r.contact?r.weight:0),0),precision.effectiveN/2);Object.assign(precision,{effectiveN:precision.effectiveN/2,ci,halfWidth:(ci[1]-ci[0])/2,met:(ci[1]-ci[0])/2<=precision.target,method:'Conservative pair-count Wilson interval'});}
   const cumulative={},traceValues=(r,key)=>key==='probability'?Array.from(r.hazardSeries,v=>coagProbability(v,p.alpha??1)):Array.from(key==='screened'?r.screenedSeries:key==='hazard'?r.hazardSeries:r.componentSeries[key]);
   for(const key of cumulativeKeys){
    const c=Object.fromEntries(['mean','median','sd','se','p05','p25','p75','p95','ciLow','ciHigh'].map(k=>[k,[]])),histories=runs.map(r=>traceValues(r,key));
    for(let b=0;b<241;b++){const q=p.modelVersion===2?weightedStats(histories.map(v=>v[b]),runs.map(r=>r.weight)):stats(histories.map(v=>v[b]));for(const k of ['mean','median','sd','se','p05','p25','p75','p95'])c[k].push(q[k]);c.ciLow.push(q.ci?.[0]??null);c.ciHigh.push(q.ci?.[1]??null);if(b%60===0)yield {completed:(species*3+o+1)*n,total:6*n,phase:'statistics'};}
    c.traces=histories.slice(0,12);c.convergence=[...new Set([Math.max(1,Math.floor(n/4)),Math.max(1,Math.floor(n/2)),n])].map(count=>({count,...summary(runs.slice(0,count),r=>r[key])}));cumulative[key]=c;
   }
   const fates=Object.fromEntries(['wall','outlet','inlet','censored','numerical'].map(k=>{const subset=runs.filter(r=>r.fate===k),count=subset.length,fraction=p.modelVersion===2?subset.reduce((s,r)=>s+r.weight,0):count/n,neff=(p.modelVersion===2?1/runs.reduce((s,r)=>s+r.weight*r.weight,0):n)/(paired?2:1);return [k,{count,fraction,ci:wilson(fraction*neff,neff)}];}));
   rows.push({species,d,offset,n,residence,precision,branchFractions:fractions,branchCounts:counts,paths:runs.filter(r=>r.points).map(r=>({points:r.points,branch:r.branchIndex,fate:r.fate})),samples,stationShear,stationCoagulation,distribution,cumulative,fates,time:summary(runs,r=>r.t),half:summary(runs,r=>r.half),quarter:summary(runs,r=>r.quarter),nh:summary(runs,r=>turnover(r.half,.5,Gref,field.peakNu,d)),nq:summary(runs,r=>turnover(r.quarter,.25,Gref,field.peakNu,d)),growth:runs.some(r=>r.growth)?{numberLoss:summary(runs,r=>r.growth.numberLoss),meanDiameter:summary(runs,r=>r.growth.meanDiameter),maxMassError:Math.max(...runs.map(r=>r.growth.massError)),maxOverflow:Math.max(...runs.map(r=>r.growth.overflowMassFraction))}:null,wall:runs.reduce((s,r)=>s+(r.contact?r.weight:0),0),wallCI:precision.ci,outlet:fates.outlet.fraction,inlet:fates.inlet.fraction,censored:fates.censored.fraction,maxRe:Math.max(...runs.map(r=>r.maxRe)),steps:Math.max(...runs.map(r=>r.steps)),curve:Array.from(mean),sd,pathFate:first.fate,path:first.points,meanPath:p.noise?track(field,p,d,offset,Gref,1,false,true).points:first.points});
  }
  return {ok:true,field,Gref,rows,etaHalf:Math.sqrt(field.peakNu/(Gref/2)),etaQuarter:Math.sqrt(field.peakNu/(Gref/4)),pressure:f,settings:{...p}};
 }
 function simulate(g,p,f,options={}){const it=simulation(g,p,f,options);let step;do{step=it.next();}while(!step.done);return step.value;}
 async function simulateAsync(g,p,f,{cancelled=()=>false,onProgress=()=>{},field,groups}={}){const it=simulation(g,p,f,{field,groups});while(true){if(cancelled()){it.return();return {ok:false,cancelled:true};}const step=it.next();if(step.done)return step.value;onProgress(step.value);await new Promise(resolve=>{if(typeof MessageChannel==='undefined'){setTimeout(resolve,0);return;}const c=new MessageChannel();c.port1.onmessage=()=>{c.port1.close();c.port2.close();resolve();};c.port2.postMessage(0);});}}
 return {numerics,dispersion,population,evolvePopulation,particle,nozzle,flow,geometry,makeField,stationSections,sampleStationCrossings,track,crossing,turnover,wilson,simulate,simulateAsync,stats,weightedStats,quantile,collisionKernel,collisionState,collisions,coagProbability,integrateRate};
})();
if(typeof module!=='undefined')module.exports=Exposure;

'use strict';
// Execution-only acceleration. Cache keys include the build, complete geometry,
// every setting, and the solved flow. Never round coordinates or physics inputs.
const ExposureRuntime=(()=>{
 function key(value){if(Array.isArray(value))return '['+value.map(key).join(',')+']';if(value&&typeof value==='object')return '{'+Object.keys(value).sort().map(k=>JSON.stringify(k)+':'+key(value[k])).join(',')+'}';return value===undefined?'undefined':Object.is(value,-0)?'-0':JSON.stringify(value);}
 class LRU{
  constructor(limit=3,bytes=48*1024*1024){this.limit=limit;this.budget=bytes;this.bytes=0;this.map=new Map();}
  get(k){const entry=this.map.get(k);if(!entry)return;this.map.delete(k);this.map.set(k,entry);return entry.value;}
  set(k,value,bytes=1){if(this.map.has(k)){this.bytes-=this.map.get(k).bytes;this.map.delete(k);}if(bytes>this.budget)return false;while(this.map.size&&(this.map.size>=this.limit||this.bytes+bytes>this.budget)){const k=this.map.keys().next().value;this.bytes-=this.map.get(k).bytes;this.map.delete(k);}this.map.set(k,{value,bytes});this.bytes+=bytes;return true;}
  clear(){this.map.clear();this.bytes=0;}
 }
 function estimate(result){let bytes=0;for(const row of result.rows){const {samples,...rest}=row;bytes+=JSON.stringify(rest).length*4;const sizes=[samples[0],samples.at(-1),samples[Math.floor(samples.length/2)]].filter(Boolean).map(s=>JSON.stringify(s).length);bytes+=Math.max(0,...sizes)*samples.length*4;}return bytes;}
 function create(options={}){
  const E=options.model||(typeof Exposure!=='undefined'?Exposure:null),W=options.WorkerClass??globalThis.Worker,build=options.buildHash||'unversioned',url=options.workerURL||'/model-assets/exposure-worker.js?v='+build;
  const cores=Math.max(1,Math.floor(options.hardwareConcurrency??globalThis.navigator?.hardwareConcurrency??2)),memory=options.deviceMemory??globalThis.navigator?.deviceMemory??8,maxWorkers=Math.max(1,Math.min(4,cores-1||1,memory<=4?2:4,options.maxWorkers??4));
  const cache=new LRU(options.cacheEntries??3,options.cacheBytes??(memory<=4?48:96)*1024*1024),prepared=new LRU(2),slots=[];let active=null,sequence=0,workersEnabled=!!W,closed=false;
  const counters={runs:0,cacheHits:0,workersStarted:0,cancelled:0,preparedHits:0};
  function remove(slot){clearTimeout(slot.timer);slot.worker.terminate();const i=slots.indexOf(slot);if(i>=0)slots.splice(i,1);}
  function cancel(){sequence++;if(!active)return;const a=active;active=null;counters.cancelled++;for(const slot of [...slots])if(slot.job?.run===a)remove(slot);a.resolve({ok:false,cancelled:true});}
  function fail(run,error,infrastructure=false){if(active!==run)return;active=null;for(const slot of [...slots])if(slot.job?.run===run)remove(slot);error.infrastructure=infrastructure;run.reject(error);}
  function dispatch(slot){const job=slot.job;if(!job||!slot.ready||job.sent)return;job.sent=true;const {run,group}=job,message={id:run.id,group,g:run.g,p:run.p,f:run.f,key:run.key,build};
   if(run.wallData&&!slot.cad.has(run.g.key)){message.wallData=run.wallData;slot.cad.add(run.g.key);}try{slot.worker.postMessage(message);}catch(error){fail(run,error,true);}
  }
  function assign(slot,run){if(active!==run||!run.queue.length)return;slot.job={run,group:run.queue.shift(),sent:false};dispatch(slot);}
  function spawn(){const slot={worker:new W(url),ready:false,job:null,cad:new Set()};slots.push(slot);counters.workersStarted++;
   slot.timer=setTimeout(()=>{if(slot.job)fail(slot.job.run,Error('Parallel engine did not start'),true);else remove(slot);},15000);
   slot.worker.onerror=e=>{e.preventDefault?.();if(slot.job)fail(slot.job.run,Error(e.message||'Parallel engine unavailable'),true);else remove(slot);};
   slot.worker.onmessageerror=()=>{if(slot.job)fail(slot.job.run,Error('Parallel result could not be transferred'),true);};
   slot.worker.onmessage=event=>{const m=event.data;if(m.ready){clearTimeout(slot.timer);if(m.build!==build){if(slot.job)fail(slot.job.run,Error('Model version changed; reload for parallel processing'),true);else remove(slot);return;}slot.ready=true;dispatch(slot);return;}
    const job=slot.job,run=job?.run;if(!run||active!==run||m.id!==run.id||m.group!==job.group)return;
    if(m.error){fail(run,Error(m.error));return;}
    if(m.progress){run.counts[m.group]=m.progress.completed;run.onProgress({completed:run.counts.reduce((a,b)=>a+b,0),total:6*run.limit,phase:m.progress.phase,workers:run.workers});return;}
    if(!m.result?.ok){fail(run,Error(m.result?.reason||'Particle model did not return a valid result'));return;}
    run.rows[m.group]=m.result.rows[0];run.meta=m.result;run.counts[m.group]=m.result.rows[0].n;slot.job=null;run.done++;
    if(run.done===6){active=null;const {rows,meta}=run;run.resolve({ok:true,rows,Gref:meta.Gref,etaHalf:meta.etaHalf,etaQuarter:meta.etaQuarter,workers:run.workers});}
    else assign(slot,run);
   };return slot;
  }
  async function simulate(g,p,f,{onProgress=()=>{},cancelled=()=>false}={}){
   if(closed)throw Error('Calculation engine was disposed');cancel();const started=performance.now(),id=++sequence,k=key({build,g,p,f});counters.runs++;
   let field=prepared.get(k);if(field)counters.preparedHits++;else{field=E.makeField(g,p,f);prepared.set(k,field);}
   const hit=cache.get(k);if(hit){counters.cacheHits++;return {...structuredClone(hit),field,pressure:f,settings:{...p},execution:{mode:'cache',workers:0,cacheHit:true,milliseconds:performance.now()-started}};}
   let result,fallbackReason=null;
   if(workersEnabled){try{const count=Math.min(maxWorkers,p.wallModel==='cad'?2:4);while(slots.length>count)remove(slots.at(-1));result=await new Promise((resolve,reject)=>{const run={id,key:k,g,p,f,wallData:options.wallData?.(g.key),queue:[0,1,2,3,4,5],rows:Array(6),counts:Array(6).fill(0),limit:p.noise?(p.sampling==='precision'?(p.maxHistories??4096):p.n):(field.geo.branches?.length||1),done:0,workers:count,resolve,reject,onProgress};active=run;try{while(slots.length<count)spawn();for(const slot of slots)assign(slot,run);}catch(error){fail(run,error,true);}});}catch(error){if(!error.infrastructure)throw error;workersEnabled=false;fallbackReason=error.message;}}
   if(result?.cancelled||cancelled()||id!==sequence)return {ok:false,cancelled:true};
   if(!result){result=await E.simulateAsync(g,p,f,{field,cancelled:()=>cancelled()||id!==sequence,onProgress});if(result.cancelled||cancelled()||id!==sequence)return {ok:false,cancelled:true};}
   if(!result.ok)return result;const data={ok:true,rows:result.rows,Gref:result.Gref,etaHalf:result.etaHalf,etaQuarter:result.etaQuarter};const bytes=estimate(data);if(bytes<=cache.budget)cache.set(k,structuredClone(data),bytes);
   return {...data,field,pressure:f,settings:{...p},execution:{mode:result.workers?'parallel':'single',workers:result.workers||1,cacheHit:false,milliseconds:performance.now()-started,...(fallbackReason?{fallbackReason}:{})}};
  }
  function dispose(){cancel();closed=true;for(const slot of [...slots])remove(slot);cache.clear();prepared.clear();}
  return {simulate,cancel,dispose,clearCache:()=>cache.clear(),stats:()=>({...counters,workers:slots.length,maxWorkers,cachedRuns:cache.map.size,cacheBytes:cache.bytes})};
 }
 return {key,LRU,estimate,create};
})();
if(typeof module!=='undefined')module.exports=ExposureRuntime;

'use strict';
// A group retains its original sequence, seed, branch weighting, adaptive
// stopping and summation order. Groups are merged by index, never finish order.
const workerFields=new ExposureRuntime.LRU(2);
self.onmessage=async event=>{
 const {id,group,g,p,f,key,wallData,build}=event.data;
 try{
  if(build!==EDUCTOR_BUILD.codeSHA256)throw Error('Worker / page model version mismatch');
  if(wallData)EductorWalls.register(g.key,wallData);
  let field=workerFields.get(key);if(!field){field=Exposure.makeField(g,p,f);workerFields.set(key,field);}
  const limit=p.noise?(p.sampling==='precision'?(p.maxHistories??4096):p.n):(field.geo.branches?.length||1);
  const r=await Exposure.simulateAsync(g,p,f,{field,groups:[group],onProgress:q=>self.postMessage({id,group,progress:{completed:q.phase==='statistics'?q.completed/(group+1):q.completed-group*limit,phase:q.phase}})});
  self.postMessage({id,group,result:{ok:r.ok,reason:r.reason,rows:r.rows,Gref:r.Gref,etaHalf:r.etaHalf,etaQuarter:r.etaQuarter}});
 }catch(error){self.postMessage({id,group,error:error.message});}
};
self.postMessage({ready:true,build:EDUCTOR_BUILD.codeSHA256});
