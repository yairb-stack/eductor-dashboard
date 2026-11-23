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
