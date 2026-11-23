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
